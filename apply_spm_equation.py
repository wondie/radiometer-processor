"""
Apply a hyperspectral band-ratio SPM (Suspended Particulate Matter) algorithm to a
multi-band remote sensing reflectance (Rrs) GeoTIFF and write the result as a new
GeoTIFF.

    X = ((Rrs550 / Rrs631) ** 3.7412904 / (Rrs550 / Rrs717) ** 3.2499739) ** 0.3693965
        * ((Rrs717 / Rrs713) ** 4.6759044 / (Rrs704 / Rrs631) ** 3.727386)

    SPM (mg/L) = 33.869 * X + 16.895

The source cube may carry more or fewer bands than the literature the algorithm was
published against (e.g. 272 bands here vs. 270 in some source papers — the Headwall
Nano-Hyperspec collection this data comes from used two different sensor
configurations, so band count varies by acquisition), so bands are never addressed by
literature index. Instead each band's wavelength is read from its GDAL
description/metadata (band names such as "HSIRrs550", "704.32 nm", or an ENVI
"wavelength" list are all supported), or — if the file carries no such metadata, as is
the case for the plain GeoTIFF exports in this project — from a linear 400-1000 nm
model matching the Nano-Hyperspec's native spectral range. The closest matching band
in *this* file is used either way, with the match reported (and flagged if it falls in
the sensor's noisier first/last ~50 nm) so the choice can be audited.

Only the five required bands are ever read from disk (not the full 272-band cube),
and the raster is processed row-stripe by row-stripe sized to the source's native
block height, so memory use stays flat regardless of image size.

Usage:
    python apply_spm_equation.py D:/HSIRrs/2018-03.tif
    python apply_spm_equation.py D:/HSIRrs/2018-03.tif -o E:/out/2018-03_SPM.tif
    python apply_spm_equation.py D:/HSIRrs/2018-03.tif --list-bands
    python apply_spm_equation.py D:/HSIRrs/2018-03.tif --band-map "550=203,631=229,717=253,713=252,704=246"

    # Folder mode: process every .tif in D:/HSIRrs, skipping ones that already
    # have output; each file gets its own optimal-drive choice unless
    # --output-dir is given. Any file's failure (e.g. no wavelength match) is
    # logged and skipped rather than stopping the rest of the batch.
    python apply_spm_equation.py D:/HSIRrs
    python apply_spm_equation.py D:/HSIRrs --output-dir E:/HSIRrs_output
    python apply_spm_equation.py D:/HSIRrs --force   # reprocess even if output exists

Requires the GDAL Python bindings (`osgeo.gdal`), e.g. via OSGeo4W, conda-forge
("conda install -c conda-forge gdal"), or a matching GDAL wheel for this Python/OS.
"""

import argparse
import os
import re
import subprocess
import sys
from typing import Dict, List, Optional

import numpy as np
from osgeo import gdal

gdal.UseExceptions()

# ---------------------------------------------------------------------------
# Algorithm
# ---------------------------------------------------------------------------
REQUIRED_WAVELENGTHS_NM = (550.0, 631.0, 717.0, 713.0, 704.0)
WAVELENGTH_MATCH_TOLERANCE_NM = 4.0  # literature (270 bands) vs. file (272 bands) drift

# Fallback band-center model, used only when the file carries no wavelength metadata
# at all (confirmed to be the case for the 2018-03.tif cube: every band description is
# empty, no ENVI domain, no sidecar .hdr). Per the official Headwall Nano-Hyperspec
# spec for this collection, the sensor natively records exactly 270 bands spanning
# 400-1000 nm (~2.2 nm/band) — this matches the 270-band literature exactly. This
# particular cube has 272 bands, i.e. 2 more than the sensor produces, and band 1 is
# empirically confirmed flat all-zero (STATISTICS_MIN=MAX=MEAN=0), so it is treated as
# a leading non-spectral/dark band. Bands 2-271 are assumed to be the sensor's 270 real
# spectral bands spanning 400-1000 nm; band 272 is left unmapped (its purpose is
# unknown — a trailing dummy band or something else — rather than guessed at).
DEFAULT_LINEAR_START_NM = 400.0
DEFAULT_LINEAR_END_NM = 1000.0
DEFAULT_DUMMY_BANDS = 1
DEFAULT_SPECTRAL_BAND_COUNT = 270

# Per the sensor spec, the first/last ~50 nm of the 400-1000 nm range are less
# reliable (higher noise, reduced accuracy). Bands matched inside this margin are
# flagged, not rejected, since it's just a data-quality caveat.
RELIABLE_RANGE_NM = (450.0, 950.0)


def linear_wavelengths(
    n_bands: int,
    start_nm: float,
    end_nm: float,
    dummy_bands: int,
    spectral_band_count: Optional[int] = None,
) -> Dict[int, float]:
    """Evenly-spaced band-center wavelengths for a run of `spectral_band_count` bands
    starting right after `dummy_bands` leading non-spectral bands. Any bands beyond
    dummy_bands + spectral_band_count (e.g. an unidentified trailing band) are left
    unmapped rather than guessed at. Defaults spectral_band_count to "all remaining
    bands" if not given.
    """
    if spectral_band_count is None:
        spectral_band_count = n_bands - dummy_bands
    if spectral_band_count < 2:
        return {}
    step = (end_nm - start_nm) / (spectral_band_count - 1)
    return {dummy_bands + i + 1: start_nm + i * step for i in range(spectral_band_count)}


def compute_spm(rrs: Dict[float, np.ndarray]) -> "tuple[np.ndarray, np.ndarray]":
    r550, r631, r717, r713, r704 = (rrs[w] for w in REQUIRED_WAVELENGTHS_NM)
    with np.errstate(divide="ignore", invalid="ignore"):
        term_a = (r550 / r631) ** 3.7412904 / (r550 / r717) ** 3.2499739
        term_b = term_a ** 0.3693965
        term_c = (r717 / r713) ** 4.6759044 / (r704 / r631) ** 3.727386
        x = term_b * term_c
        spm = 33.869 * x + 16.895
    return x.astype(np.float32), spm.astype(np.float32)


# ---------------------------------------------------------------------------
# Band wavelength discovery
# ---------------------------------------------------------------------------
_NM_SUFFIX_PATTERN = re.compile(r"([-+]?\d*\.\d+|\d+)\s*(?:nm|nanometers?)", re.IGNORECASE)
_DECIMAL_PATTERN = re.compile(r"([-+]?\d+\.\d+)")
_PLAUSIBLE_WL_PATTERN = re.compile(r"(\d{3,4})(?:\.\d+)?")


def _extract_wavelength(text: str) -> Optional[float]:
    """Pull a wavelength (nm) out of a band description/name.

    Handles, in priority order: an explicit "<value> nm" suffix, a bare decimal
    number (wavelengths carry decimals, band indices usually don't), and finally
    an embedded 3-4 digit integer in a plausible optical range (handles names
    like "HSIRrs550" with no separator or unit at all).
    """
    if not text:
        return None
    match = _NM_SUFFIX_PATTERN.search(text)
    if match:
        return float(match.group(1))
    match = _DECIMAL_PATTERN.search(text)
    if match:
        return float(match.group(1))
    matches = list(_PLAUSIBLE_WL_PATTERN.finditer(text))
    if matches:
        value = float(matches[-1].group(1))
        if 350.0 <= value <= 2500.0:
            return value
    return None


def get_band_wavelengths(ds: gdal.Dataset) -> Dict[int, float]:
    """Return {1-based band index: wavelength_nm} using whatever the file exposes."""
    wavelengths: Dict[int, float] = {}
    for i in range(1, ds.RasterCount + 1):
        band = ds.GetRasterBand(i)
        meta = band.GetMetadata()
        candidate = meta.get("wavelength") or meta.get("Wavelength") or meta.get("WAVELENGTH")
        wl = _extract_wavelength(candidate) if candidate else None
        if wl is None:
            wl = _extract_wavelength(band.GetDescription())
        if wl is not None:
            wavelengths[i] = wl

    if not wavelengths:
        # Fall back to a dataset-level ENVI-style "wavelength = {v1, v2, ...}" list.
        for domain in ("ENVI", ""):
            dmeta = ds.GetMetadata(domain)
            wl_str = dmeta.get("wavelength") or dmeta.get("Wavelength")
            if wl_str:
                values = [float(v) for v in re.findall(r"[-+]?\d*\.?\d+", wl_str)]
                if len(values) == ds.RasterCount:
                    wavelengths = {i + 1: v for i, v in enumerate(values)}
                    break

    if wavelengths and max(wavelengths.values()) < 20.0:
        # Values look like micrometers (e.g. 0.55) rather than nanometers.
        wavelengths = {i: v * 1000.0 for i, v in wavelengths.items()}

    return wavelengths


def match_bands(
    wavelengths: Dict[int, float], targets: "tuple[float, ...]", tolerance: float
) -> Dict[float, int]:
    matched: Dict[float, int] = {}
    for target in targets:
        best_band, best_diff = None, None
        for band, wl in wavelengths.items():
            diff = abs(wl - target)
            if best_diff is None or diff < best_diff:
                best_band, best_diff = band, diff
        if best_band is None or best_diff > tolerance:
            raise ValueError(
                f"No band within {tolerance} nm of {target} nm "
                f"(closest: {best_diff} nm). Use --band-map to override."
            )
        matched[target] = best_band
    return matched


def load_wavelengths_file(path: str) -> Dict[int, float]:
    """Read one wavelength (nm) per line/row, in 1-based band order.

    Accepts a plain text/CSV file such as the band-center calibration list that
    ships with most pushbroom hyperspectral cameras (Resonon, Headwall, Specim,
    ...). Each non-empty line may be just a number, or "index,wavelength" /
    "index wavelength"; the last numeric token on the line is used.
    """
    wavelengths: Dict[int, float] = {}
    with open(path, "r", encoding="utf-8-sig") as fh:
        band_idx = 0
        for line in fh:
            tokens = re.findall(r"[-+]?\d*\.?\d+", line.strip())
            if not tokens:
                continue  # blank line or header row (e.g. "Wavelength,FWHM")
            band_idx += 1
            wavelengths[band_idx] = float(tokens[-1])
    return wavelengths


# Real per-band wavelength calibration for Headwall Nano-Hyperspec unit "nhs-008",
# extracted from D:/HSIRrs/MS.x839.000.00L.hdr's ENVI "wavelength" list (272 bands,
# 397.118-998.738 nm, ~2.22 nm spacing). This is measured calibration data, not a
# model, so it's preferred over the linear fallback whenever the open file's band
# count matches — but only then, since it's specific to this one sensor unit.
BUNDLED_CALIBRATION_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "nhs008_wavelengths.txt")


def load_bundled_calibration(n_bands: int) -> Dict[int, float]:
    if not os.path.exists(BUNDLED_CALIBRATION_PATH):
        return {}
    wavelengths = load_wavelengths_file(BUNDLED_CALIBRATION_PATH)
    return wavelengths if len(wavelengths) == n_bands else {}


def parse_band_map(band_map_str: str) -> Dict[float, int]:
    mapping: Dict[float, int] = {}
    for pair in band_map_str.split(","):
        wl_str, idx_str = pair.split("=")
        mapping[float(wl_str.strip())] = int(idx_str.strip())
    missing = set(REQUIRED_WAVELENGTHS_NM) - set(mapping)
    if missing:
        raise ValueError(f"--band-map is missing wavelengths: {sorted(missing)}")
    return mapping


def print_band_matches(matches: Dict[float, int], wavelengths: Dict[int, float]) -> None:
    print("Band matches (target nm -> file band # @ actual nm):")
    for target in REQUIRED_WAVELENGTHS_NM:
        band_idx = matches[target]
        actual = wavelengths.get(band_idx, target)
        flag = ""
        if not (RELIABLE_RANGE_NM[0] <= actual <= RELIABLE_RANGE_NM[1]):
            flag = f"  [WARNING: outside sensor's reliable {RELIABLE_RANGE_NM[0]:.0f}-{RELIABLE_RANGE_NM[1]:.0f} nm range]"
        print(f"  {target:>6.1f} nm -> band {band_idx:>4d} @ {actual:>8.2f} nm (delta {actual - target:+.2f} nm){flag}")


# ---------------------------------------------------------------------------
# Optimal output location (SSD/NVMe aware)
# ---------------------------------------------------------------------------
def _drive_letter(path: str) -> str:
    return os.path.splitdrive(os.path.abspath(path))[0].upper()


_ssd_cache: Dict[str, Optional[bool]] = {}


def _is_ssd(drive_letter: str) -> Optional[bool]:
    """Best-effort SSD/NVMe detection on Windows. Returns True/False, or None if
    it can't be determined (non-Windows, no permissions, network/virtual drive).
    Memoized per drive letter so batch runs over many files on the same drive
    don't spawn a PowerShell process per file."""
    if drive_letter in _ssd_cache:
        return _ssd_cache[drive_letter]
    if os.name != "nt":
        _ssd_cache[drive_letter] = None
        return None
    ps_cmd = (
        "$ErrorActionPreference='Stop'; "
        f"$vol = Get-Volume -DriveLetter '{drive_letter.rstrip(':')}'; "
        "$part = Get-Partition -Volume $vol; "
        "(Get-PhysicalDisk -DeviceNumber $part.DiskNumber).MediaType"
    )
    try:
        result = subprocess.run(
            ["powershell", "-NoProfile", "-NonInteractive", "-Command", ps_cmd],
            capture_output=True, text=True, timeout=10,
        )
        media = result.stdout.strip()
        is_ssd = True if media == "SSD" else False if media == "HDD" else None
    except Exception:
        is_ssd = None
    _ssd_cache[drive_letter] = is_ssd
    return is_ssd


def choose_output_path(input_path: str, explicit_output: Optional[str], relative_dir: str = "") -> str:
    """`relative_dir`, if given, is the input file's subfolder path relative to a
    batch root — used only by the alternate-drive branch below, which otherwise
    pools every file from a recursive batch into one shared folder and would
    collide on same-named files from different subfolders (e.g. two different
    flights both containing a "Flight1_FL1_rf").
    """
    if explicit_output:
        os.makedirs(os.path.dirname(os.path.abspath(explicit_output)) or ".", exist_ok=True)
        return explicit_output

    input_drive = _drive_letter(input_path)
    ssd = _is_ssd(input_drive)
    base = os.path.splitext(os.path.basename(input_path))[0]

    if ssd is False:
        # Spinning disk: reading input and writing output on the same physical
        # disk causes seek contention between the two streams. Prefer a
        # different drive if one is available.
        for candidate in ("C:", "D:", "E:", "F:"):
            if candidate != input_drive and os.path.isdir(candidate + "\\"):
                out_dir = os.path.join(candidate + "\\", "HSIRs_output", relative_dir)
                os.makedirs(out_dir, exist_ok=True)
                print(f"{input_drive} looks like an HDD; writing output to {candidate} to avoid seek contention.")
                return os.path.join(out_dir, f"{base}_SPM.tif")

    # SSD/NVMe, or no alternative drive available, or detection was inconclusive:
    # solid-state media has no seek penalty, so keeping input and output on the
    # same drive is simplest and avoids a cross-drive copy.
    out_dir = os.path.join(os.path.dirname(os.path.abspath(input_path)), "output")
    os.makedirs(out_dir, exist_ok=True)
    if ssd:
        print(f"{input_drive} detected as SSD/NVMe; writing output alongside the input on {input_drive}.")
    else:
        print(f"Could not determine media type for {input_drive}; defaulting to writing alongside the input.")
    return os.path.join(out_dir, f"{base}_SPM.tif")


# ---------------------------------------------------------------------------
# Raster processing
# ---------------------------------------------------------------------------
def read_scaled_band(band: gdal.Band, x0: int, y0: int, xsize: int, ysize: int) -> np.ndarray:
    arr = band.ReadAsArray(x0, y0, xsize, ysize).astype(np.float64)
    scale = band.GetScale()
    offset = band.GetOffset()
    if scale not in (None, 1.0):
        arr *= scale
    if offset not in (None, 0.0):
        arr += offset
    return arr


def create_output_dataset(driver: gdal.Driver, path: str, x_size: int, y_size: int, bands: int):
    base_options = ["TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256", "BIGTIFF=IF_SAFER", "NUM_THREADS=ALL_CPUS"]
    for compress_options in (
        ["COMPRESS=ZSTD", "ZSTD_LEVEL=9", "PREDICTOR=3"],
        ["COMPRESS=DEFLATE", "ZLEVEL=6", "PREDICTOR=3"],
        [],
    ):
        try:
            return driver.Create(path, x_size, y_size, bands, gdal.GDT_Float32, base_options + compress_options)
        except RuntimeError:
            continue
    raise RuntimeError(f"Could not create output dataset at {path}")


def resolve_band_matches(
    ds: gdal.Dataset,
    tolerance: float,
    band_map_str: Optional[str] = None,
    wavelengths_file: Optional[str] = None,
    use_linear_fallback: bool = True,
    linear_start_nm: float = DEFAULT_LINEAR_START_NM,
    linear_end_nm: float = DEFAULT_LINEAR_END_NM,
    dummy_bands: int = DEFAULT_DUMMY_BANDS,
    spectral_band_count: Optional[int] = DEFAULT_SPECTRAL_BAND_COUNT,
) -> "tuple[Dict[float, int], Dict[int, float]]":
    """Resolve REQUIRED_WAVELENGTHS_NM to 1-based band indices in `ds`, trying (in
    order): --band-map, --wavelengths-file, embedded file metadata, the bundled
    nhs-008 calibration (if band count matches), then the linear model fallback.
    Shared by apply_spm_equation.py and extract_spm_bands.py so both always agree
    on which physical bands the algorithm's wavelengths resolve to.
    """
    n_bands = ds.RasterCount
    if band_map_str:
        matches = parse_band_map(band_map_str)
        wavelengths = {idx: wl for wl, idx in matches.items()}
        return matches, wavelengths

    if wavelengths_file:
        wavelengths = load_wavelengths_file(wavelengths_file)
    else:
        wavelengths = get_band_wavelengths(ds)
        if not wavelengths:
            wavelengths = load_bundled_calibration(n_bands)
            if wavelengths:
                print(
                    f"No wavelength metadata on any band; using the bundled nhs-008 "
                    f"sensor calibration ({n_bands} real measured wavelengths from "
                    f"MS.x839.000.00L.hdr) since the band count matches."
                )
    if not wavelengths and n_bands == len(REQUIRED_WAVELENGTHS_NM):
        # A file with exactly as many bands as the algorithm needs is almost
        # certainly extract_spm_bands.py's own output (possibly re-mosaicked
        # afterward, which typically drops the per-band SetDescription() labels
        # that would otherwise let get_band_wavelengths() find these directly).
        # Trust the canonical band order that script always writes in, rather
        # than running the 270-band linear model against a 5-band file (which
        # produces nonsense band indices like "band 69" and crashes).
        wavelengths = {i + 1: wl for i, wl in enumerate(REQUIRED_WAVELENGTHS_NM)}
        order_str = ", ".join(f"{i + 1}=Rrs{wl:g}" for i, wl in enumerate(REQUIRED_WAVELENGTHS_NM))
        print(
            f"No wavelength metadata found, but this file has exactly "
            f"{n_bands} bands, matching the number of Rrs bands the algorithm "
            f"needs. Assuming it's already extract_spm_bands.py's output in its "
            f"canonical band order: {order_str}. Use --band-map to override if "
            f"this file's bands are ordered differently."
        )
    if not wavelengths and use_linear_fallback:
        wavelengths = linear_wavelengths(n_bands, linear_start_nm, linear_end_nm, dummy_bands, spectral_band_count)
        n_spectral = spectral_band_count if spectral_band_count is not None else n_bands - dummy_bands
        first_band, last_band = dummy_bands + 1, dummy_bands + n_spectral
        step = (linear_end_nm - linear_start_nm) / (n_spectral - 1)
        extra = f"; bands {last_band + 1}-{n_bands} left unmapped (unidentified)" if last_band < n_bands else ""
        print(
            f"No wavelength metadata found on any band (no descriptions, no ENVI "
            f"domain, no sidecar .hdr). Falling back to the Headwall Nano-Hyperspec's "
            f"native 270-band linear model: bands {first_band}-{last_band} spanning "
            f"{linear_start_nm:.1f}-{linear_end_nm:.1f} nm ({step:.3f} nm/band){extra}. "
            f"If the original ENVI .hdr for this flight line still exists (before "
            f"GeoTIFF conversion/mosaicking), point --wavelengths-file at its "
            f"wavelength list instead for an exact mapping. Override with "
            f"--start-nm/--end-nm/--dummy-bands/--spectral-bands, --wavelengths-file, "
            f"or --band-map."
        )
    if not wavelengths:
        raise RuntimeError(
            "No wavelength metadata found on any band, and no --wavelengths-file/"
            "--band-map was given, and --no-linear-fallback was set. Supply one of:\n"
            "  --wavelengths-file path/to/band_centers.txt   "
            "(one wavelength in nm per line, 1-based band order)\n"
            '  --band-map "550=<idx>,631=<idx>,717=<idx>,713=<idx>,704=<idx>"  '
            "(exact 1-based band indices, if already known)"
        )
    matches = match_bands(wavelengths, REQUIRED_WAVELENGTHS_NM, tolerance)
    return matches, wavelengths


def process(
    input_path: str,
    output_path: str,
    tolerance: float,
    band_map_str: Optional[str] = None,
    wavelengths_file: Optional[str] = None,
    use_linear_fallback: bool = True,
    linear_start_nm: float = DEFAULT_LINEAR_START_NM,
    linear_end_nm: float = DEFAULT_LINEAR_END_NM,
    dummy_bands: int = DEFAULT_DUMMY_BANDS,
    spectral_band_count: Optional[int] = DEFAULT_SPECTRAL_BAND_COUNT,
) -> None:
    ds = gdal.Open(input_path, gdal.GA_ReadOnly)
    x_size, y_size, n_bands = ds.RasterXSize, ds.RasterYSize, ds.RasterCount
    print(f"Opened {input_path}: {x_size} x {y_size}, {n_bands} bands")

    matches, wavelengths = resolve_band_matches(
        ds, tolerance, band_map_str, wavelengths_file, use_linear_fallback,
        linear_start_nm, linear_end_nm, dummy_bands, spectral_band_count,
    )
    print_band_matches(matches, wavelengths)

    src_bands = [ds.GetRasterBand(matches[wl]) for wl in REQUIRED_WAVELENGTHS_NM]
    nodata = src_bands[0].GetNoDataValue()

    _, block_y = src_bands[0].GetBlockSize()
    stripe_height = max(block_y, 256)  # sequential, cache-friendly stripes

    driver = gdal.GetDriverByName("GTiff")
    out_ds = create_output_dataset(driver, output_path, x_size, y_size, 1)
    out_ds.SetGeoTransform(ds.GetGeoTransform())
    out_ds.SetProjection(ds.GetProjection())

    out_spm_band = out_ds.GetRasterBand(1)
    out_spm_band.SetDescription("SPM (mg/L)")
    out_spm_band.SetNoDataValue(float("nan"))

    for y0 in range(0, y_size, stripe_height):
        rows = min(stripe_height, y_size - y0)
        arrays = [read_scaled_band(b, 0, y0, x_size, rows) for b in src_bands]
        rrs = dict(zip(REQUIRED_WAVELENGTHS_NM, arrays))
        _, spm = compute_spm(rrs)

        if nodata is not None:
            invalid = np.zeros((rows, x_size), dtype=bool)
            for arr in arrays:
                invalid |= arr == nodata
            spm[invalid] = np.nan

        out_spm_band.WriteArray(spm, 0, y0)

    out_spm_band.FlushCache()
    out_ds.FlushCache()
    print(f"Wrote output to {output_path}")


# ---------------------------------------------------------------------------
# Folder-mode raster discovery
# ---------------------------------------------------------------------------
# Folder names we generate outputs into ourselves; never descend into or list
# files from these when scanning for inputs.
OUTPUT_FOLDER_NAMES = {"output", "spm_bands", "hsirs_output"}

# Known Headwall/hyperspectral flight-line naming conventions, adapted from the
# QGIS batchmosaicker plugin's hyperspectral_detector.py (built and tested against
# this same HSIRrs collection, which mixes several naming schemes across flights).
# These are extension-less ENVI binary rasters, so a sibling ".hdr" is required
# too — the pattern match alone just avoids false-positives on unrelated
# extension-less files GDAL happens to be able to open.
_HSI_FILENAME_PATTERNS = [re.compile(p, re.IGNORECASE) for p in (
    r".*_rf$",                        # *_rf (reflectance-corrected)
    r"^FL\d+$",                       # FL1, FL2, ...
    r"^Flight\d+[a-zA-Z]?_FL\d+$",    # Flight1_FL1
    r"^Flight\d+[a-zA-Z]?_F\d+$",     # Flight10_F1
    r"^(North|South|East|West)_FL$",
    r".*Full_Flight.*",
    r".*flig[ht]+Path\d+.*",          # flightPath / fligthPath (typo in some data)
    r"^Box\d+_FL\d+$",
    r"^Box\d+_\d+_FL\d+$",
    r"^\d+[A-Z]?_FL\d+$",
    r"^F\d+_FL\d+$",
)]


def _is_envi_raster_without_extension(path: str) -> bool:
    """True for Headwall's per-flight-line ENVI exports: an extension-less raw
    binary raster (e.g. "Flight1_FL1_rf") with a sibling "<name>.hdr" header
    (required for GDAL's ENVI driver to open it) whose filename also matches a
    known flight-line naming convention.
    """
    filename = os.path.basename(path)
    if os.path.splitext(filename)[1]:
        return False  # has some extension; handled by the .tif branch or excluded
    if not os.path.exists(path + ".hdr"):
        return False
    return any(pattern.match(filename) for pattern in _HSI_FILENAME_PATTERNS)


def _prefer_rf_over_raw(paths: List[str]) -> List[str]:
    """Drop raw (radiance/DN) files that have a "<name>_rf" reflectance-corrected
    sibling in the same folder — the SPM equation needs Rrs, not raw radiance, so
    processing the raw version alongside its _rf counterpart would silently waste
    time on data the algorithm can't use correctly.
    """
    by_dir: Dict[str, set] = {}
    for p in paths:
        by_dir.setdefault(os.path.dirname(p), set()).add(os.path.basename(p))
    return [
        p for p in paths
        if os.path.basename(p).endswith("_rf")
        or f"{os.path.basename(p)}_rf" not in by_dir.get(os.path.dirname(p), ())
    ]


def find_raster_inputs(root_dir: str, recursive: bool = True) -> List[str]:
    """Find hyperspectral raster inputs under root_dir: GeoTIFFs (excluding this
    project's own *_SPM.tif/*_bands.tif outputs) and extension-less ENVI rasters
    that match a known flight-line naming convention and have a sibling .hdr (see
    _is_envi_raster_without_extension). Prefers _rf files over their raw
    counterpart (see _prefer_rf_over_raw). Skips descending into this project's
    own output folders. Returns full paths, sorted.
    """
    def is_wanted(filename: str, full_path: str) -> bool:
        lower = filename.lower()
        if lower.endswith(".tif"):
            return not (lower.endswith("_spm.tif") or lower.endswith("_bands.tif"))
        return _is_envi_raster_without_extension(full_path)

    results: List[str] = []
    if recursive:
        for dirpath, dirnames, filenames in os.walk(root_dir):
            dirnames[:] = [d for d in dirnames if d.lower() not in OUTPUT_FOLDER_NAMES]
            for filename in filenames:
                full_path = os.path.join(dirpath, filename)
                if is_wanted(filename, full_path):
                    results.append(full_path)
    else:
        for filename in os.listdir(root_dir):
            full_path = os.path.join(root_dir, filename)
            if os.path.isfile(full_path) and is_wanted(filename, full_path):
                results.append(full_path)
    return sorted(_prefer_rf_over_raw(results))


def batch_process(
    input_dir: str,
    output_dir: Optional[str],
    tolerance: float,
    band_map_str: Optional[str],
    wavelengths_file: Optional[str],
    use_linear_fallback: bool,
    linear_start_nm: float,
    linear_end_nm: float,
    dummy_bands: int,
    spectral_band_count: Optional[int],
    force: bool,
    recursive: bool = True,
) -> None:
    """Run process() for every raster found under input_dir (see find_raster_inputs),
    skipping any that already have a corresponding output (unless force=True). One
    file's failure doesn't stop the rest.
    """
    input_paths = find_raster_inputs(input_dir, recursive)
    if not input_paths:
        print(f"No raster files found under {input_dir}")
        return

    processed, skipped, failed = [], [], []
    for input_path in input_paths:
        base = os.path.splitext(os.path.basename(input_path))[0]
        relative_dir = os.path.relpath(os.path.dirname(input_path), input_dir)
        relative_dir = "" if relative_dir == "." else relative_dir

        if output_dir:
            # Mirror the input's subfolder structure under output_dir — flattening
            # would collide, since different flights reuse identical filenames
            # (e.g. "Flight1_FL1_rf" appears under multiple date folders).
            nested_dir = os.path.join(output_dir, relative_dir) if relative_dir else output_dir
            os.makedirs(nested_dir, exist_ok=True)
            output_path = os.path.join(nested_dir, f"{base}_SPM.tif")
        else:
            output_path = choose_output_path(input_path, None, relative_dir=relative_dir)

        if os.path.exists(output_path) and not force:
            print(f"[skip] {input_path}: output already exists at {output_path}")
            skipped.append(input_path)
            continue

        print(f"[run]  {input_path} -> {output_path}")
        try:
            process(
                input_path, output_path, tolerance, band_map_str, wavelengths_file,
                use_linear_fallback, linear_start_nm, linear_end_nm, dummy_bands,
                spectral_band_count,
            )
            processed.append(input_path)
        except Exception as exc:
            print(f"[FAIL] {input_path}: {exc}")
            failed.append(input_path)

    print(
        f"\nBatch complete: {len(processed)} processed, {len(skipped)} skipped "
        f"(output already existed), {len(failed)} failed."
    )
    if failed:
        print("Failed files:", ", ".join(failed))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", nargs="?", default="D:/HSIRrs/2018-03.tif", help="A single hyperspectral raster (GeoTIFF or ENVI), or a folder to process every raster found in it (recursively by default).")
    parser.add_argument("-o", "--output", default=None, help="Explicit output path for a single input file; if omitted, chosen automatically.")
    parser.add_argument("--output-dir", default=None, help="Directory for all outputs when input is a folder. If omitted, the optimal drive is chosen per file automatically.")
    parser.add_argument("--force", action="store_true", help="Reprocess files even if a corresponding output already exists (folder mode only).")
    parser.add_argument("--no-recursive", action="store_true", help="Only scan the top level of the input folder instead of descending into subfolders (folder mode only).")
    parser.add_argument("--tolerance", type=float, default=WAVELENGTH_MATCH_TOLERANCE_NM, help="Max wavelength mismatch (nm) allowed when matching bands.")
    parser.add_argument("--band-map", default=None, help='Manual override, e.g. "550=203,631=229,717=253,713=252,704=246".')
    parser.add_argument("--wavelengths-file", default=None, help="Text/CSV file with one band-center wavelength (nm) per line, in 1-based band order.")
    parser.add_argument("--start-nm", type=float, default=DEFAULT_LINEAR_START_NM, help="Linear-fallback first spectral band's wavelength (nm), used only when no metadata/--wavelengths-file/--band-map is available.")
    parser.add_argument("--end-nm", type=float, default=DEFAULT_LINEAR_END_NM, help="Linear-fallback last spectral band's wavelength (nm).")
    parser.add_argument("--dummy-bands", type=int, default=DEFAULT_DUMMY_BANDS, help="Number of leading non-spectral bands to skip in the linear fallback (default: 1, matching the confirmed all-zero band 1).")
    parser.add_argument("--spectral-bands", type=int, default=DEFAULT_SPECTRAL_BAND_COUNT, help="Number of real spectral bands in the linear fallback, starting right after --dummy-bands (default: 270, the Nano-Hyperspec's documented native band count). Remaining bands, if any, are left unmapped.")
    parser.add_argument("--no-linear-fallback", action="store_true", help="Disable the linear 400-1000nm Headwall Nano-Hyperspec assumption; fail instead if no wavelength metadata is found.")
    parser.add_argument("--list-bands", action="store_true", help="Print discovered band wavelengths and matches, then exit without writing output.")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        sys.exit(f"Input path not found: {args.input}")

    gdal.SetCacheMax(1024 * 1024 * 1024)  # 1 GiB block cache

    if os.path.isdir(args.input):
        if args.list_bands:
            sys.exit("--list-bands works on a single file, not a folder.")
        output_dir = args.output_dir or args.output
        if args.output_dir and args.output and args.output_dir != args.output:
            print(f"Both -o/--output and --output-dir given for a folder input; using --output-dir ({args.output_dir}).")
        batch_process(
            args.input, output_dir, args.tolerance, args.band_map, args.wavelengths_file,
            use_linear_fallback=not args.no_linear_fallback,
            linear_start_nm=args.start_nm, linear_end_nm=args.end_nm, dummy_bands=args.dummy_bands,
            spectral_band_count=args.spectral_bands, force=args.force, recursive=not args.no_recursive,
        )
        return

    if args.list_bands:
        ds = gdal.Open(args.input, gdal.GA_ReadOnly)
        wavelengths = load_wavelengths_file(args.wavelengths_file) if args.wavelengths_file else get_band_wavelengths(ds)
        if not wavelengths and not args.wavelengths_file:
            wavelengths = load_bundled_calibration(ds.RasterCount)
            if wavelengths:
                print(f"Using bundled nhs-008 sensor calibration ({ds.RasterCount} real measured wavelengths)")
        if not wavelengths and not args.no_linear_fallback:
            wavelengths = linear_wavelengths(ds.RasterCount, args.start_nm, args.end_nm, args.dummy_bands, args.spectral_bands)
            last_band = args.dummy_bands + args.spectral_bands
            print(f"Using linear fallback: bands {args.dummy_bands + 1}-{last_band} span {args.start_nm:.1f}-{args.end_nm:.1f} nm")
        if not wavelengths:
            sys.exit("No wavelength metadata found on any band. Try --wavelengths-file.")
        for band_idx in sorted(wavelengths):
            print(f"  band {band_idx:>4d}: {wavelengths[band_idx]:>8.2f} nm")
        matches = match_bands(wavelengths, REQUIRED_WAVELENGTHS_NM, args.tolerance)
        print_band_matches(matches, wavelengths)
        return

    output_path = choose_output_path(args.input, args.output)
    process(
        args.input, output_path, args.tolerance, args.band_map, args.wavelengths_file,
        use_linear_fallback=not args.no_linear_fallback,
        linear_start_nm=args.start_nm, linear_end_nm=args.end_nm, dummy_bands=args.dummy_bands,
        spectral_band_count=args.spectral_bands,
    )


if __name__ == "__main__":
    main()
