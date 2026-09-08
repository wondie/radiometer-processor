"""
Compute each mission's mean SPM (mg/L) from its raster (D:/HSIRrs/output/
{year}_{month:02d}_SPM.tif, produced by apply_spm_equation.py), and write that
mission-wide mean into the "UAS SPM" column of every sites_SPM.xlsx row whose
date falls in that mission — as a new workbook, so the original is untouched.

A mission-wide mean (rather than sampling the pixel under each site's exact
coordinate) sidesteps a real problem with per-point sampling: many sample-site
coordinates land just outside the UAS flight swath (the Headwall processing
notes for this data explicitly flag residual GPS/IMU georeferencing uncertainty
with no tie points to correct it), so a meaningful fraction of sites would get
NoData if sampled individually. The per-mission mean is also exactly what
plot/box_plot_line_hyper.py's hardcoded `stats[mission]['mean']` values already
use — this script recomputes them fresh from the current rasters and cross-checks
against those hardcoded values so drift between the two is caught, not silent.

Rows with no matching raster for their mission, or with a 2021 date, are dropped
from the output workbook entirely rather than left with a stale/blank UAS SPM.

Usage:
    python sample_spm_at_sites.py
    python sample_spm_at_sites.py --spm-dir D:/HSIRrs/output --xlsx path/to/sites_SPM.xlsx
"""

import argparse
import datetime
import os
import sys
from typing import Dict, Optional

import numpy as np
from osgeo import gdal

gdal.UseExceptions()

DEFAULT_SPM_DIR = "D:/HSIRrs/output"
DEFAULT_XLSX = "G:/Other computers/My Laptop/dissertation/SPM_Hyperspectral/data/sites_SPM.xlsx"
DEFAULT_OUTPUT_XLSX = "G:/Other computers/My Laptop/dissertation/SPM_Hyperspectral/data/sites_SPM_hyper.xlsx"
DEFAULT_SHEET_NAME = "sites_SPM_time_series"
DATE_COLUMN_NAME = "Date"
UAS_COLUMN_NAME = "UAS SPM"

# Reference means already computed and hardcoded in plot/box_plot_line_hyper.py's
# `stats` dict, used only to sanity-check the freshly computed values below.
REFERENCE_MEANS = {
    "2018-03": 48.49413536058343,
    "2018-05": 50.66261584940906,
    "2018-06": 53.106575937526294,
    "2018-07": 48.4147312334088,
    "2018-12": 54.76202082100937,
    "2019-06": 39.861843893526085,
    "2019-07": 37.61548303090227,
}


# A handful of pixels in some missions (e.g. 2018-07, 2018-12) carry physically
# impossible values up to ~1e26 mg/L — near-zero-denominator blowups in
# compute_spm()'s unclipped power/ratio terms, affecting well under 0.2% of
# pixels but large enough to destroy a simple mean. Excluded here rather than
# silently baked into the source rasters, since fixing that upstream means
# deciding how compute_spm() should handle it (NaN out, clip, ...) and
# regenerating every raster — a separate decision from computing this mean.
MAX_PLAUSIBLE_SPM = 1000.0


def mission_mean(raster_path: str) -> "tuple[float, int, int]":
    """Mean of valid, physically plausible pixels (0 <= SPM <= MAX_PLAUSIBLE_SPM),
    computed in row-stripes so memory use stays flat regardless of raster size.
    Returns (mean, valid_pixel_count, excluded_outlier_count).
    """
    ds = gdal.Open(raster_path, gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    x_size, y_size = ds.RasterXSize, ds.RasterYSize
    _, block_y = band.GetBlockSize()
    stripe_height = max(block_y, 512)

    total, count, excluded = 0.0, 0, 0
    for y0 in range(0, y_size, stripe_height):
        rows = min(stripe_height, y_size - y0)
        arr = band.ReadAsArray(0, y0, x_size, rows)
        not_nan = arr[~np.isnan(arr)]
        plausible = (not_nan >= 0) & (not_nan <= MAX_PLAUSIBLE_SPM)
        total += float(not_nan[plausible].sum())
        count += int(plausible.sum())
        excluded += not_nan.size - int(plausible.sum())

    if count == 0:
        raise RuntimeError(f"No valid pixels found in {raster_path}")
    return total / count, count, excluded


def find_spm_rasters(spm_dir: str) -> Dict[str, str]:
    """{'YYYY-MM': path} for every {year}_{month:02d}_SPM.tif in spm_dir."""
    rasters = {}
    for filename in os.listdir(spm_dir):
        if filename.endswith("_SPM.tif"):
            year_month = filename[: -len("_SPM.tif")].replace("_", "-")
            rasters[year_month] = os.path.join(spm_dir, filename)
    return rasters


# A mission's fieldwork can spill a day or two into the next calendar month
# (e.g. the July 2019 mission's last site visits landed on Aug 1-2 — confirmed
# by the workbook's own "Month-Year" column tagging those rows "2019-07-31").
# Only days this early in the month are treated as a spillover from the
# previous month's mission; anything later is a genuinely separate mission.
TRAILING_DAY_THRESHOLD = 5


def resolve_mission(date_value, available_missions) -> Optional[str]:
    year_month = date_value.strftime("%Y-%m")
    if year_month in available_missions:
        return year_month
    if date_value.day <= TRAILING_DAY_THRESHOLD:
        previous_month = (date_value.replace(day=1) - datetime.timedelta(days=1)).strftime("%Y-%m")
        if previous_month in available_missions:
            return previous_month
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--spm-dir", default=DEFAULT_SPM_DIR, help="Folder containing {year}_{month:02d}_SPM.tif rasters.")
    parser.add_argument("--xlsx", default=DEFAULT_XLSX, help="Existing sites_SPM.xlsx to read (never modified in place).")
    parser.add_argument("--output-xlsx", default=DEFAULT_OUTPUT_XLSX, help="Path for the new workbook with updated UAS SPM values.")
    parser.add_argument("--sheet", default=DEFAULT_SHEET_NAME, help="Worksheet name in the xlsx.")
    parser.add_argument("--in-place", action="store_true", help="Overwrite --xlsx directly instead of writing --output-xlsx.")
    args = parser.parse_args()

    import openpyxl

    if not os.path.isdir(args.spm_dir):
        sys.exit(f"SPM raster folder not found: {args.spm_dir}")
    if not os.path.exists(args.xlsx):
        sys.exit(f"Workbook not found: {args.xlsx}")

    rasters = find_spm_rasters(args.spm_dir)
    if not rasters:
        sys.exit(f"No *_SPM.tif rasters found in {args.spm_dir}")

    print(f"Computing mission means from {len(rasters)} raster(s)...")
    means: Dict[str, float] = {}
    for year_month, path in sorted(rasters.items()):
        mean, valid_count, excluded_count = mission_mean(path)
        means[year_month] = mean
        reference = REFERENCE_MEANS.get(year_month)
        flag = ""
        if reference is not None and abs(mean - reference) > 0.01:
            flag = f"  [DIFFERS from box_plot_line_hyper.py's hardcoded {reference:.4f} by {mean - reference:+.4f}]"
        outlier_note = f", excluded {excluded_count} outlier pixel(s) (>{'{:.0f}'.format(MAX_PLAUSIBLE_SPM)} mg/L)" if excluded_count else ""
        print(f"  {year_month}: {mean:.4f} mg/L over {valid_count} valid pixels ({os.path.basename(path)}{outlier_note}){flag}")

    wb = openpyxl.load_workbook(args.xlsx)
    ws = wb[args.sheet]

    header = [cell.value for cell in next(ws.iter_rows(min_row=1, max_row=1))]
    try:
        date_col = header.index(DATE_COLUMN_NAME) + 1
        uas_col = header.index(UAS_COLUMN_NAME) + 1
    except ValueError as exc:
        sys.exit(f"Expected columns {DATE_COLUMN_NAME!r} and {UAS_COLUMN_NAME!r} not found in header {header}: {exc}")

    updated_counts: Dict[str, int] = {}
    removed_counts: Dict[str, int] = {}
    rows_to_delete = []
    for row_idx in range(2, ws.max_row + 1):
        date_value = ws.cell(row=row_idx, column=date_col).value
        is_2021 = date_value is not None and date_value.year == 2021
        mission = resolve_mission(date_value, means) if date_value is not None else None
        mean = means.get(mission) if mission else None

        if mean is None or is_2021:
            reason = (date_value.strftime("%Y-%m") if date_value is not None else "no date")
            removed_counts[reason] = removed_counts.get(reason, 0) + 1
            rows_to_delete.append(row_idx)
            continue

        ws.cell(row=row_idx, column=uas_col, value=mean)
        updated_counts[mission] = updated_counts.get(mission, 0) + 1

    # Delete bottom-to-top so earlier indices in rows_to_delete stay valid.
    for row_idx in reversed(rows_to_delete):
        ws.delete_rows(row_idx, 1)

    output_path = args.xlsx if args.in_place else args.output_xlsx
    wb.save(output_path)

    print(f"\nUpdated {sum(updated_counts.values())} row(s) across {len(updated_counts)} mission(s) in {output_path}")
    for year_month, count in sorted(updated_counts.items()):
        print(f"  {year_month}: {count} row(s) -> {means[year_month]:.4f} mg/L")
    if removed_counts:
        print(f"\nRemoved {sum(removed_counts.values())} row(s) with no matching raster or a 2021 date:")
        for reason, count in sorted(removed_counts.items()):
            print(f"  {reason}: {count} row(s)")


if __name__ == "__main__":
    main()
