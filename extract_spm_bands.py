"""
Extract just the 5 Rrs bands that apply_spm_equation.py's algorithm uses
(550, 631, 717, 713, 704 nm) from a hyperspectral raster — GeoTIFF or ENVI — or
every raster found in a folder (recursively by default), and write them as a
small multi-band GeoTIFF into a "spm_bands" folder.

Band resolution (embedded metadata -> bundled nhs-008 calibration -> linear
fallback -> --wavelengths-file/--band-map overrides) and raster discovery (*.tif
plus extension-less ENVI files with a sibling .hdr, e.g. Headwall's per-flight-line
exports) are shared with apply_spm_equation.py, so both scripts always agree on
which physical bands the algorithm's wavelengths point to and which files count
as inputs.

Usage:
    python extract_spm_bands.py D:/HSIRrs/2018-03.tif
    python extract_spm_bands.py D:/HSIRrs                     # every raster, recursive
    python extract_spm_bands.py D:/HSIRrs --force
    python extract_spm_bands.py D:/HSIRrs/2018-03.tif -o D:/HSIRrs/spm_bands/2018-03_bands.tif
    python extract_spm_bands.py D:/HSIRrs/MS.x839.000.0018-data -o D:/HSI/2018_06   # folder input: -o = --output-dir
"""

import argparse
import os
import sys
from typing import Optional

import numpy as np
from osgeo import gdal

from apply_spm_equation import (
    DEFAULT_DUMMY_BANDS,
    DEFAULT_LINEAR_END_NM,
    DEFAULT_LINEAR_START_NM,
    DEFAULT_SPECTRAL_BAND_COUNT,
    REQUIRED_WAVELENGTHS_NM,
    WAVELENGTH_MATCH_TOLERANCE_NM,
    create_output_dataset,
    find_raster_inputs,
    print_band_matches,
    read_scaled_band,
    resolve_band_matches,
)

gdal.UseExceptions()

SPM_BANDS_FOLDER_NAME = "spm_bands"


def spm_bands_output_path(input_path: str, output_dir: Optional[str], relative_dir: str = "") -> str:
    """`relative_dir`, if given, is the input file's subfolder path relative to a
    batch root, mirrored under output_dir so that a recursive folder run doesn't
    flatten same-named files from different subfolders into one location.
    """
    base = os.path.splitext(os.path.basename(input_path))[0]
    if output_dir is None:
        out_dir = os.path.join(os.path.dirname(os.path.abspath(input_path)), SPM_BANDS_FOLDER_NAME)
    else:
        out_dir = os.path.join(output_dir, relative_dir) if relative_dir else output_dir
    os.makedirs(out_dir, exist_ok=True)
    return os.path.join(out_dir, f"{base}_bands.tif")


def extract_bands(
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
    out_ds = create_output_dataset(driver, output_path, x_size, y_size, len(REQUIRED_WAVELENGTHS_NM))
    out_ds.SetGeoTransform(ds.GetGeoTransform())
    out_ds.SetProjection(ds.GetProjection())

    out_bands = []
    for i, wl in enumerate(REQUIRED_WAVELENGTHS_NM, start=1):
        band = out_ds.GetRasterBand(i)
        src_idx = matches[wl]
        band.SetDescription(f"Rrs{wl:g} (source band {src_idx} @ {wavelengths[src_idx]:.2f} nm)")
        band.SetNoDataValue(float("nan"))
        out_bands.append(band)

    for y0 in range(0, y_size, stripe_height):
        rows = min(stripe_height, y_size - y0)
        arrays = [read_scaled_band(b, 0, y0, x_size, rows) for b in src_bands]

        if nodata is not None:
            invalid = np.zeros((rows, x_size), dtype=bool)
            for arr in arrays:
                invalid |= arr == nodata
            for arr in arrays:
                arr[invalid] = np.nan

        for out_band, arr in zip(out_bands, arrays):
            out_band.WriteArray(arr.astype(np.float32), 0, y0)

    for out_band in out_bands:
        out_band.FlushCache()
    out_ds.FlushCache()
    print(f"Wrote output to {output_path}")


def batch_extract(
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
    """Run extract_bands() for every raster found under input_dir (see
    find_raster_inputs), skipping any that already have a corresponding output
    (unless force=True). One file's failure doesn't stop the rest.
    """
    input_paths = find_raster_inputs(input_dir, recursive)
    if not input_paths:
        print(f"No raster files found under {input_dir}")
        return

    processed, skipped, failed = [], [], []
    for input_path in input_paths:
        relative_dir = os.path.relpath(os.path.dirname(input_path), input_dir)
        relative_dir = "" if relative_dir == "." else relative_dir
        output_path = spm_bands_output_path(input_path, output_dir, relative_dir=relative_dir)

        if os.path.exists(output_path) and not force:
            print(f"[skip] {input_path}: output already exists at {output_path}")
            skipped.append(input_path)
            continue

        print(f"[run]  {input_path} -> {output_path}")
        try:
            extract_bands(
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
    parser.add_argument("-o", "--output", default=None, help="Explicit output path for a single input file; if omitted, written to a 'spm_bands' subfolder next to the input.")
    parser.add_argument("--output-dir", default=None, help="Folder for all outputs (single file or folder mode). If omitted, uses '<input dir>/spm_bands'.")
    parser.add_argument("--force", action="store_true", help="Reprocess files even if a corresponding output already exists (folder mode only).")
    parser.add_argument("--no-recursive", action="store_true", help="Only scan the top level of the input folder instead of descending into subfolders (folder mode only).")
    parser.add_argument("--tolerance", type=float, default=WAVELENGTH_MATCH_TOLERANCE_NM, help="Max wavelength mismatch (nm) allowed when matching bands.")
    parser.add_argument("--band-map", default=None, help='Manual override, e.g. "550=203,631=229,717=253,713=252,704=246".')
    parser.add_argument("--wavelengths-file", default=None, help="Text/CSV file with one band-center wavelength (nm) per line, in 1-based band order.")
    parser.add_argument("--start-nm", type=float, default=DEFAULT_LINEAR_START_NM, help="Linear-fallback first spectral band's wavelength (nm).")
    parser.add_argument("--end-nm", type=float, default=DEFAULT_LINEAR_END_NM, help="Linear-fallback last spectral band's wavelength (nm).")
    parser.add_argument("--dummy-bands", type=int, default=DEFAULT_DUMMY_BANDS, help="Number of leading non-spectral bands to skip in the linear fallback.")
    parser.add_argument("--spectral-bands", type=int, default=DEFAULT_SPECTRAL_BAND_COUNT, help="Number of real spectral bands in the linear fallback, starting right after --dummy-bands.")
    parser.add_argument("--no-linear-fallback", action="store_true", help="Disable the linear 400-1000nm Headwall Nano-Hyperspec assumption; fail instead if no wavelength metadata is found.")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        sys.exit(f"Input path not found: {args.input}")

    gdal.SetCacheMax(1024 * 1024 * 1024)  # 1 GiB block cache

    if os.path.isdir(args.input):
        output_dir = args.output_dir or args.output
        if args.output_dir and args.output and args.output_dir != args.output:
            print(f"Both -o/--output and --output-dir given for a folder input; using --output-dir ({args.output_dir}).")
        batch_extract(
            args.input, output_dir, args.tolerance, args.band_map, args.wavelengths_file,
            not args.no_linear_fallback, args.start_nm, args.end_nm, args.dummy_bands,
            args.spectral_bands, args.force, recursive=not args.no_recursive,
        )
        return

    if args.output:
        output_path = args.output
        os.makedirs(os.path.dirname(os.path.abspath(output_path)) or ".", exist_ok=True)
    else:
        output_path = spm_bands_output_path(args.input, args.output_dir)

    extract_bands(
        args.input, output_path, args.tolerance, args.band_map, args.wavelengths_file,
        not args.no_linear_fallback, args.start_nm, args.end_nm, args.dummy_bands,
        args.spectral_bands,
    )


if __name__ == "__main__":
    main()
