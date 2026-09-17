"""
Apply a MicaSense-specific band-ratio SPM (Suspended Particulate Matter) algorithm to
per-band MicaSense RedEdge orthomosaics and write the result as a single-band GeoTIFF.

    x = (R_rs717 / R_rs668) ** 3.2363 + (R_rs842 * R_rs668 ** 4.2625) / (R_rs717 ** 5.2195)
    SPM (mg/L) = 46.56 * x - 32.63

Unlike apply_spm_equation.py (one hyperspectral cube with bands identified by
wavelength), MicaSense RedEdge missions in D:/Micasense/Mission are exported as one
single-band GeoTIFF per camera band, named "<mission>_<band>.tif" where <band> is the
camera's fixed 1-5 firmware band number (1=Blue, 2=Green, 3=Red ~668nm, 4=NIR ~842nm,
5=RedEdge ~717nm — see MicaSense's IMG_XXXX_B<N> file-naming documentation). So bands
here are identified by that fixed index, not by reading wavelength metadata.

Per-band files for the same mission are not always pixel-identical in extent (e.g. for
2018_06 the NIR file is 73 columns narrower than the Red/RedEdge files despite sharing
the same origin and pixel size), so this script intersects the three bands' georeferenced
extents and reads only the overlapping window from each, rather than assuming identical
raster dimensions.

Usage:
    python apply_spm_equation_micasense.py D:/Micasense/Mission
    python apply_spm_equation_micasense.py D:/Micasense/Mission --force
    python apply_spm_equation_micasense.py D:/Micasense/Mission --output-dir D:/Micasense/Mission/output
    python apply_spm_equation_micasense.py D:/Micasense/Mission/2018_06_3.tif
    python apply_spm_equation_micasense.py D:/Micasense/Mission/2018_06_3.tif -o D:/out/2018_06_SPM.tif

Requires the GDAL Python bindings (`osgeo.gdal`).
"""

import argparse
import os
import re
import sys
from typing import Dict, List, Optional

import numpy as np
from osgeo import gdal

from apply_spm_equation import choose_output_path, create_output_dataset, read_scaled_band

gdal.UseExceptions()

# Fixed MicaSense RedEdge firmware band numbers (not detected from metadata — see
# module docstring). Only the three bands the algorithm needs are opened.
BAND_NAMES = {1: "Blue", 2: "Green", 3: "Red", 4: "NIR", 5: "RedEdge"}
RED_BAND, NIR_BAND, REDEDGE_BAND = 3, 4, 5
REQUIRED_BANDS = (RED_BAND, NIR_BAND, REDEDGE_BAND)

_BAND_FILENAME_PATTERN = re.compile(r"^(?P<mission>.+)_(?P<band>[1-5])\.tif$", re.IGNORECASE)


def parse_band_filename(path: str) -> Optional["tuple[str, int]"]:
    """Return (mission_name, band_number) for a "<mission>_<1-5>.tif" filename, or
    None if it doesn't match that convention."""
    match = _BAND_FILENAME_PATTERN.match(os.path.basename(path))
    if not match:
        return None
    return match.group("mission"), int(match.group("band"))


def find_missions(root_dir: str, recursive: bool = True) -> Dict[str, Dict[int, str]]:
    """Group "<mission>_<band>.tif" files under root_dir by mission. Returns
    {mission_path_prefix: {band_number: file_path}}, where mission_path_prefix is
    the full path with the "_<band>.tif" suffix stripped.
    """
    missions: Dict[str, Dict[int, str]] = {}

    def scan(dirpath: str, filenames: List[str]) -> None:
        for filename in filenames:
            parsed = parse_band_filename(filename)
            if parsed is None:
                continue
            mission, band = parsed
            key = os.path.join(dirpath, mission)
            missions.setdefault(key, {})[band] = os.path.join(dirpath, filename)

    if recursive:
        for dirpath, _dirnames, filenames in os.walk(root_dir):
            scan(dirpath, filenames)
    else:
        scan(root_dir, os.listdir(root_dir))

    return missions


def compute_spm(r717: np.ndarray, r668: np.ndarray, r842: np.ndarray) -> "tuple[np.ndarray, np.ndarray]":
    with np.errstate(divide="ignore", invalid="ignore"):
        x = (r717 / r668) ** 3.2363 + (r842 * r668 ** 4.2625) / (r717 ** 5.2195)
        spm = 46.56 * x - 32.63
    return x.astype(np.float32), spm.astype(np.float32)


def common_grid(datasets: Dict[int, "gdal.Dataset"]) -> "tuple[tuple, int, int, Dict[int, tuple]]":
    """Compute the georeferenced intersection of same-projection rasters that share
    (approximately) the same origin and pixel size but may differ slightly in extent.
    Returns (output_geotransform, width, height, {band: (xoff, yoff)}), where each
    offset is in that band's own dataset's pixel space.
    """
    ref_gt = next(iter(datasets.values())).GetGeoTransform()
    pixel_w, pixel_h = ref_gt[1], ref_gt[5]  # pixel_h is negative (north-up)

    lefts, rights, tops, bottoms = [], [], [], []
    for ds in datasets.values():
        gt = ds.GetGeoTransform()
        left, top = gt[0], gt[3]
        right = left + ds.RasterXSize * gt[1]
        bottom = top + ds.RasterYSize * gt[5]
        lefts.append(left)
        rights.append(right)
        tops.append(top)
        bottoms.append(bottom)

    left_common, right_common = max(lefts), min(rights)
    top_common, bottom_common = min(tops), max(bottoms)
    width = round((right_common - left_common) / pixel_w)
    height = round((top_common - bottom_common) / -pixel_h)
    if width <= 0 or height <= 0:
        raise RuntimeError("Bands do not spatially overlap.")

    offsets: Dict[int, tuple] = {}
    for band, ds in datasets.items():
        gt = ds.GetGeoTransform()
        xoff = round((left_common - gt[0]) / gt[1])
        yoff = round((top_common - gt[3]) / gt[5])
        offsets[band] = (xoff, yoff)

    out_gt = (left_common, pixel_w, 0.0, top_common, 0.0, pixel_h)
    return out_gt, width, height, offsets


def process_mission(band_paths: Dict[int, str], output_path: str) -> None:
    datasets = {band: gdal.Open(band_paths[band], gdal.GA_ReadOnly) for band in REQUIRED_BANDS}
    for band, ds in datasets.items():
        print(f"  {BAND_NAMES[band]} (band {band}): {band_paths[band]} — {ds.RasterXSize} x {ds.RasterYSize}")

    gt, x_size, y_size, offsets = common_grid(datasets)
    print(f"  Common overlap: {x_size} x {y_size}")

    src_bands = {band: ds.GetRasterBand(1) for band, ds in datasets.items()}
    nodata = {band: b.GetNoDataValue() for band, b in src_bands.items()}

    _, block_y = src_bands[RED_BAND].GetBlockSize()
    stripe_height = max(block_y, 256)

    driver = gdal.GetDriverByName("GTiff")
    out_ds = create_output_dataset(driver, output_path, x_size, y_size, 1)
    out_ds.SetGeoTransform(gt)
    out_ds.SetProjection(datasets[RED_BAND].GetProjection())

    out_band = out_ds.GetRasterBand(1)
    out_band.SetDescription("SPM (mg/L)")
    out_band.SetNoDataValue(float("nan"))

    for y0 in range(0, y_size, stripe_height):
        rows = min(stripe_height, y_size - y0)
        arrays = {}
        for band in REQUIRED_BANDS:
            xoff, yoff = offsets[band]
            arrays[band] = read_scaled_band(src_bands[band], xoff, yoff + y0, x_size, rows)

        invalid = np.zeros((rows, x_size), dtype=bool)
        for band, arr in arrays.items():
            if nodata[band] is not None:
                invalid |= arr == nodata[band]

        _, spm = compute_spm(arrays[REDEDGE_BAND], arrays[RED_BAND], arrays[NIR_BAND])
        invalid |= ~np.isfinite(spm)
        spm[invalid] = np.nan

        out_band.WriteArray(spm, 0, y0)

    out_band.FlushCache()
    out_ds.FlushCache()
    print(f"  Wrote output to {output_path}")


def batch_process(input_dir: str, output_dir: Optional[str], force: bool, recursive: bool = True) -> None:
    missions = find_missions(input_dir, recursive)
    if not missions:
        print(f"No '<mission>_<1-5>.tif' files found under {input_dir}")
        return

    processed, skipped, failed = [], [], []
    for mission_prefix, band_paths in sorted(missions.items()):
        missing = [b for b in REQUIRED_BANDS if b not in band_paths]
        if missing:
            names = ", ".join(f"{b}={BAND_NAMES[b]}" for b in missing)
            print(f"[skip] {mission_prefix}: missing required band(s) {names}")
            skipped.append(mission_prefix)
            continue

        base = os.path.basename(mission_prefix)
        if output_dir:
            relative_dir = os.path.relpath(os.path.dirname(mission_prefix), input_dir)
            relative_dir = "" if relative_dir == "." else relative_dir
            nested_dir = os.path.join(output_dir, relative_dir) if relative_dir else output_dir
            os.makedirs(nested_dir, exist_ok=True)
            output_path = os.path.join(nested_dir, f"{base}_SPM.tif")
        else:
            output_path = choose_output_path(mission_prefix, None)

        if os.path.exists(output_path) and not force:
            print(f"[skip] {mission_prefix}: output already exists at {output_path}")
            skipped.append(mission_prefix)
            continue

        print(f"[run]  {mission_prefix} -> {output_path}")
        try:
            process_mission(band_paths, output_path)
            processed.append(mission_prefix)
        except Exception as exc:
            print(f"[FAIL] {mission_prefix}: {exc}")
            failed.append(mission_prefix)

    print(
        f"\nBatch complete: {len(processed)} processed, {len(skipped)} skipped "
        f"(output already existed or incomplete), {len(failed)} failed."
    )
    if failed:
        print("Failed missions:", ", ".join(failed))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "input", nargs="?", default="D:/Micasense/Mission",
        help="A folder containing MicaSense per-band GeoTIFFs named '<mission>_<1-5>.tif' "
             "(processes every complete mission found, recursively by default), or a single "
             "band file (e.g. '2018_06_3.tif') to process just that one mission.",
    )
    parser.add_argument("-o", "--output", default=None, help="Explicit output path (single-mission input only).")
    parser.add_argument("--output-dir", default=None, help="Folder for all outputs (folder input only). If omitted, the optimal drive is chosen per mission automatically.")
    parser.add_argument("--force", action="store_true", help="Reprocess missions even if a corresponding output already exists (folder mode only).")
    parser.add_argument("--no-recursive", action="store_true", help="Only scan the top level of the input folder instead of descending into subfolders (folder mode only).")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        sys.exit(f"Input path not found: {args.input}")

    gdal.SetCacheMax(1024 * 1024 * 1024)  # 1 GiB block cache

    if os.path.isdir(args.input):
        batch_process(args.input, args.output_dir or args.output, args.force, recursive=not args.no_recursive)
        return

    parsed = parse_band_filename(args.input)
    if parsed is None:
        sys.exit(f"{args.input} doesn't match the '<mission>_<1-5>.tif' naming convention.")
    mission, _band = parsed
    mission_prefix = os.path.join(os.path.dirname(os.path.abspath(args.input)), mission)

    band_paths = {}
    for band in REQUIRED_BANDS:
        candidate = f"{mission_prefix}_{band}.tif"
        if not os.path.exists(candidate):
            sys.exit(f"Missing required band {band} ({BAND_NAMES[band]}) file: {candidate}")
        band_paths[band] = candidate

    output_path = choose_output_path(mission_prefix, args.output)
    process_mission(band_paths, output_path)


if __name__ == "__main__":
    main()
