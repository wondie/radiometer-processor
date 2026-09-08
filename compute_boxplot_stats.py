"""
Recompute plot/box_plot_line_hyper.py's hardcoded `stats` dict from the actual
SPM rasters in D:/HSIRrs/output, replacing the old values (which came from a
different, older raster source: E:/SPMOutput/SPM_{month}.tif, per that file's
now-dead raster_to_pixel_values()/txt_file_to_pixel_values() functions).

Uses matplotlib.cbook.boxplot_stats() on each mission's full pixel distribution
(all valid, physically-plausible pixels — see MAX_PLAUSIBLE_SPM in
sample_spm_at_sites.py for why that filter exists) so the stats are directly
comparable to whatever originally produced the hardcoded dict. `fliers` is
reset to [] before printing, matching the original dict's convention — with
hundreds of millions of pixels, the flier list itself would be enormous and
isn't meant to be embedded as a literal.

Usage:
    python compute_boxplot_stats.py
    python compute_boxplot_stats.py --spm-dir D:/HSIRrs/output
"""

import argparse
import os

import matplotlib.cbook as cbook
import numpy as np
from osgeo import gdal

from sample_spm_at_sites import MAX_PLAUSIBLE_SPM, find_spm_rasters

gdal.UseExceptions()


def load_plausible_pixels(raster_path: str) -> np.ndarray:
    """All valid, physically-plausible pixel values from band 1, as one flat array."""
    ds = gdal.Open(raster_path, gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    x_size, y_size = ds.RasterXSize, ds.RasterYSize
    _, block_y = band.GetBlockSize()
    stripe_height = max(block_y, 512)

    chunks = []
    for y0 in range(0, y_size, stripe_height):
        rows = min(stripe_height, y_size - y0)
        arr = band.ReadAsArray(0, y0, x_size, rows)
        not_nan = arr[~np.isnan(arr)]
        plausible = not_nan[(not_nan >= 0) & (not_nan <= MAX_PLAUSIBLE_SPM)]
        if plausible.size:
            chunks.append(plausible)
    return np.concatenate(chunks) if chunks else np.array([], dtype=np.float32)


def format_stats_dict(mission_stats: dict) -> str:
    entries = []
    for year_month in sorted(mission_stats.keys()):
        s = {k: (float(v) if isinstance(v, (np.floating, np.integer)) else v) for k, v in mission_stats[year_month].items()}
        entries.append(
            f"{year_month!r}: {{'label': {s['label']!r}, 'mean': {s['mean']!r},\n"
            f"                     'iqr': {s['iqr']!r}, 'cilo': {s['cilo']!r},\n"
            f"                     'cihi': {s['cihi']!r}, 'whishi': {s['whishi']!r}, 'whislo': {s['whislo']!r},\n"
            f"                     'fliers': [], 'q1': {s['q1']!r}, 'med': {s['med']!r},\n"
            f"                     'q3': {s['q3']!r}}}"
        )
    body = (",\n         ").join(entries)
    return f"stats = {{{body}\n         }}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--spm-dir", default="D:/HSIRrs/output", help="Folder containing {year}_{month:02d}_SPM.tif rasters.")
    args = parser.parse_args()

    rasters = find_spm_rasters(args.spm_dir)
    if not rasters:
        raise SystemExit(f"No *_SPM.tif rasters found in {args.spm_dir}")

    mission_stats = {}
    for year_month, path in sorted(rasters.items()):
        print(f"Loading {os.path.basename(path)}...")
        data = load_plausible_pixels(path)
        print(f"  {data.size} plausible pixels, computing boxplot stats...")
        result = cbook.boxplot_stats(data, labels=[year_month])[0]
        result["fliers"] = []

        year, month = year_month.split("-")
        result["label"] = f"{month}-01-{year}"
        mission_stats[year_month] = result
        print(f"  mean={result['mean']:.4f} med={result['med']:.4f} q1={result['q1']:.4f} q3={result['q3']:.4f} "
              f"whislo={result['whislo']:.4f} whishi={result['whishi']:.4f}")

    print("\n" + format_stats_dict(mission_stats))


if __name__ == "__main__":
    main()
