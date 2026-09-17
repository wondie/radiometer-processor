"""
Validate a published band-ratio SPM (Suspended Particulate Matter) algorithm against
independent in-situ data, and produce a combined two-panel figure:

    (a) training diagnostic - the algorithm's X term, computed on --training-path,
        plotted against that data's own measured SPM, with a locally best-fit line
        and R^2 (a sanity check on how well this algorithm's band-ratio structure
        correlates with SPM on the data at hand - not a refit of the algorithm
        itself).
    (b) validation          - measured vs. estimated SPM on --testing-path, where
        "estimated" uses the algorithm's *published* fixed slope/intercept (the
        actual equation as deployed elsewhere in this project, e.g.
        apply_spm_equation.py), with the 1:1 line and RMSE/MAE/Bias/rRMSE error
        statistics.

Panel rendering is delegated to plot/plot.py's create_scatter_plot (the same
function algorithm_generator2.py's validate_model() uses), so the figure matches
the rest of this project's SPM validation plots.

Algorithms (see ALGORITHMS below) are looked up by --algorithm:

  hyperspectral (default) - the Headwall Nano-Hyperspec equation from
      apply_spm_equation.py:
          X = ((Rrs550/Rrs631)**3.7412904 / (Rrs550/Rrs717)**3.2499739)**0.3693965
              * ((Rrs717/Rrs713)**4.6759044 / (Rrs704/Rrs631)**3.727386)
          SPM (mg/L) = 33.869 * X + 16.895

  modis - SPM (mg/L) = -17.489 * X + 61.833, where
          X = (Rrs678*Rrs412)**2.828605 / (Rrs443*Rrs443)**2.901816
              - (Rrs867/Rrs443)**3.836643
      No default training/testing data is bundled yet for this one - pass
      --training-path/--testing-path once MODIS in-situ data is available.

Usage:
    python validate_spm_algorithm.py
    python validate_spm_algorithm.py --algorithm hyperspectral \\
        --training-path data/hypers_total_rrs.csv \\
        --testing-path "G:/Other computers/My Laptop/dissertation/SPM Testing/2026/hyperspectral_testing.xlsx"
    python validate_spm_algorithm.py --algorithm modis \\
        --training-path path/to/modis_training.csv --testing-path path/to/modis_testing.xlsx
"""

import argparse
import os
import sys
from collections import OrderedDict
from typing import Dict, Optional

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold

from plot import plot

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))


def _hyperspectral_x(r: Dict[int, np.ndarray]) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        term_a = (r[550] / r[631]) ** 3.7412904 / (r[550] / r[717]) ** 3.2499739
        term_b = term_a ** 0.3693965
        term_c = (r[717] / r[713]) ** 4.6759044 / (r[704] / r[631]) ** 3.727386
        return term_b * term_c


def _modis_x(r: Dict[int, np.ndarray]) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        return (r[678] * r[412]) ** 2.828605 / (r[443] * r[443]) ** 2.901816 - (r[867] / r[443]) ** 3.836643


ALGORITHMS = {
    "hyperspectral": {
        "bands_nm": (550, 631, 717, 713, 704),
        "compute_x": _hyperspectral_x,
        "slope": 33.869,
        "intercept": 16.895,
        # add_algorithm_plot() wraps this whole string in a single pair of `$`,
        # but mathtext can't line-break inside one $...$ span - so the "$\n$"
        # here closes the first math span and reopens a second one, giving two
        # independently-parsed lines instead of one that overflows the axis.
        "equation_latex": (
            r"\left(\left(R_{rs}550/R_{rs}631\right)^{3.7412904}"
            r"/\left(R_{rs}550/R_{rs}717\right)^{3.2499739}\right)^{0.3693965}"
            "$\n$"
            r"\times\left(\left(R_{rs}717/R_{rs}713\right)^{4.6759044}"
            r"/\left(R_{rs}704/R_{rs}631\right)^{3.727386}\right)"
        ),
        # data/hypers_total_rrs.csv (this repo) is only a 30-row subset of the
        # actual 50-sample calibration set this equation's published constants
        # (33.869/16.895, R^2=0.72, F=123.3 on df=48) were fit on - using it as
        # the default understated R^2 by roughly half. This file is the real one.
        "default_training_path": r"G:\Other computers\My Laptop\dissertation\SPM_Hyperspectral\Field_Data\Hyperspectral2020\Hyperspectral_SPM_2020.csv",
        "default_testing_path": r"G:\Other computers\My Laptop\dissertation\SPM Testing\2026\hyperspectral_testing.xlsx",
        "default_testing_sheet": "Sheet2",
    },
    "modis": {
        "bands_nm": (678, 412, 443, 867),
        "compute_x": _modis_x,
        "slope": -17.489,
        "intercept": 61.833,
        "equation_latex": (
            r"\left(R_{rs}678\times R_{rs}412\right)^{2.828605}"
            r"/\left(R_{rs}443\times R_{rs}443\right)^{2.901816}"
            r"-\left(R_{rs}867/R_{rs}443\right)^{3.836643}"
        ),
        # modis.csv's "Calc" column is the algorithm's X term already computed
        # (presumably from real MODIS-band Rrs, not band-convolved hyperspectral
        # like the testing set below) - confirmed against modis-MLR_coeff.csv,
        # whose fitted (Intercept, Calc) = (61.833379, -17.489111) match the
        # published equation exactly.
        "default_training_path": r"G:\Other computers\My Laptop\dissertation\SPM_MODIS\validation\modis.csv",
        "default_training_x_column": "Calc",
        # Testing data is derived from the same hyperspectral site file used for
        # the "hyperspectral" algorithm, band-convolved down to MODIS bands via
        # the SRF below (see convolve_to_modis_bands()) rather than measured
        # directly by a MODIS sensor.
        "default_testing_path": r"G:\Other computers\My Laptop\dissertation\SPM Testing\2026\hyperspectral_testing.xlsx",
        "default_testing_sheet": "Sheet2",
        "testing_srf_path": r"G:\Other computers\My Laptop\dissertation\SPM_MODIS\Spectral_response_function_modis.csv",
    },
}

# SRF CSV column -> MODIS ocean-color band center (nm), matching this project's
# own convention from calculate_regression.R's `modis_headers` (note band12 is
# labeled 555 nm there, band16 867 nm - close to but not exactly each band's
# true nominal center - kept as-is so convolved columns line up with
# ALGORITHMS["modis"]["bands_nm"] via build_band_dict's bare-int lookup).
MODIS_SRF_BAND_TO_NM = {
    "band8": 412, "band9": 443, "band10": 488, "band11": 531,
    "band12": 555, "band13": 667, "band14": 678, "band15": 750, "band16": 867,
}


class FixedLinearModel:
    """Adapts a published (slope, intercept) pair to the sklearn-style
    `.coef_`/`.intercept_` interface plot.add_algorithm_plot() reads when
    labeling panel (a) - without pretending it was fit by sklearn.
    """

    def __init__(self, slope: float, intercept: float):
        self.coef_ = np.array([slope])
        self.intercept_ = intercept

    def predict(self, x: np.ndarray) -> np.ndarray:
        return self.coef_[0] * np.ravel(x) + self.intercept_


def build_band_dict(df: pd.DataFrame, bands_nm: "tuple[int, ...]") -> Dict[int, np.ndarray]:
    """Look up each required band as either a bare-int column (e.g. the
    hyperspectral_testing.xlsx convention) or an "Rrs<nm>" column (e.g.
    hypers_total_rrs.csv's convention).
    """
    result: Dict[int, np.ndarray] = {}
    for band in bands_nm:
        if band in df.columns:
            result[band] = df[band].to_numpy(dtype=float)
        elif f"Rrs{band}" in df.columns:
            result[band] = df[f"Rrs{band}"].to_numpy(dtype=float)
        else:
            raise KeyError(
                f"No column for band {band} nm (tried {band!r} and 'Rrs{band}'). "
                f"Available columns: {list(df.columns)[:10]}..."
            )
    return result


def _read_table(path: str, sheet_name: Optional[str] = None) -> pd.DataFrame:
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xls"):
        return pd.read_excel(path, sheet_name=sheet_name) if sheet_name else pd.read_excel(path)
    return pd.read_csv(path)


def _x_and_y_from_df(df: pd.DataFrame, bands_nm: "tuple[int, ...]", compute_x) -> "tuple[np.ndarray, np.ndarray]":
    band_dict = build_band_dict(df, bands_nm)
    x = compute_x(band_dict)
    y = df["SPM"].to_numpy(dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    return x[valid], y[valid]


def load_dataset(path: str, bands_nm: "tuple[int, ...]", compute_x, sheet_name: Optional[str] = None, x_column: Optional[str] = None) -> "tuple[np.ndarray, np.ndarray]":
    """If `x_column` is given, the file already carries the algorithm's X term
    precomputed under that column name (e.g. modis.csv's "Calc") - used as-is
    instead of recomputing X from raw Rrs band columns.
    """
    df = _read_table(path, sheet_name)
    if x_column is not None:
        x = df[x_column].to_numpy(dtype=float)
        y = df["SPM"].to_numpy(dtype=float)
        valid = np.isfinite(x) & np.isfinite(y)
        return x[valid], y[valid]
    return _x_and_y_from_df(df, bands_nm, compute_x)


def convolve_to_modis_bands(df: pd.DataFrame, srf_path: str) -> pd.DataFrame:
    """Band-convolve a per-nm hyperspectral Rrs table down to MODIS ocean-color
    bands: for each band, a weighted average of Rrs(wavelength) over that
    band's spectral response, sum(Rrs*SRF)/sum(SRF) - the SRF file's weights
    don't sum to 1 on their own, so skipping that division silently inflates
    every convolved value by ~10-16x.

    Only wavelengths present in both the SRF and `df` are used; SPM (if
    present) passes through unchanged.
    """
    srf = pd.read_csv(srf_path).set_index("Wavelength")
    wl = sorted(int(c) for c in df.columns if isinstance(c, (int, float)) and int(c) in srf.index)
    if not wl:
        raise ValueError("No wavelength columns in the site file overlap the SRF's wavelength range.")

    rrs = df[wl].to_numpy(dtype=float)
    srf_arr = srf.loc[wl].to_numpy(dtype=float)
    with np.errstate(invalid="ignore"):
        modis_rrs = (rrs @ srf_arr) / srf_arr.sum(axis=0)

    band_names = [MODIS_SRF_BAND_TO_NM[c] for c in srf.columns]
    result = pd.DataFrame(modis_rrs, columns=band_names, index=df.index)
    if "SPM" in df.columns:
        result["SPM"] = df["SPM"].to_numpy()

    # Negative Rrs is physically impossible (a raw-spectrum noise artifact, not
    # a real signal) and breaks any fractional-power term on that band, so drop
    # those sites here rather than let them silently NaN out downstream.
    negative_rows = (result[band_names] < 0).any(axis=1)
    if negative_rows.any():
        site_col = "wavelength" if "wavelength" in df.columns else None
        names = df.loc[negative_rows, site_col].tolist() if site_col else negative_rows[negative_rows].index.tolist()
        print(f"Dropping {negative_rows.sum()} site(s) with negative convolved Rrs (raw spectrum noise): {names}")
        result = result.loc[~negative_rows]
    return result


def build_modis_band_table(df: pd.DataFrame, srf_path: str) -> pd.DataFrame:
    """Same convolution as convolve_to_modis_bands(), but formatted for export:
    site names kept as a leading "Site" column, and MODIS band columns labeled
    "band8_412".."band16_867" (this project's convention, e.g.
    calculate_regression.R's modis_headers / WMS_Rrs_modis_SPM.xlsx) instead of
    bare nm ints, so the file is self-describing when opened on its own.
    """
    srf_cols = pd.read_csv(srf_path, nrows=0).columns.tolist()
    srf_cols = [c for c in srf_cols if c != "Wavelength"]
    modis_df = convolve_to_modis_bands(df, srf_path)
    rename = {MODIS_SRF_BAND_TO_NM[c]: f"{c}_{MODIS_SRF_BAND_TO_NM[c]}" for c in srf_cols}
    modis_df = modis_df.rename(columns=rename)
    site_col = "wavelength" if "wavelength" in df.columns else df.columns[0]
    # .loc by modis_df's (possibly convolve_to_modis_bands-filtered) index, not
    # a positional .to_numpy() straight off df - row counts can differ now that
    # negative-Rrs sites get dropped during convolution.
    modis_df.insert(0, "Site", df.loc[modis_df.index, site_col].to_numpy())
    band_cols = [rename[MODIS_SRF_BAND_TO_NM[c]] for c in srf_cols]
    ordered_cols = ["Site"] + (["SPM"] if "SPM" in modis_df.columns else []) + band_cols
    return modis_df[ordered_cols]


def flag_scale_outliers(df: pd.DataFrame, wl_range: "tuple[int, int]" = (400, 900), threshold: float = 2.0) -> None:
    """Print a warning for any row whose mean reflectance over `wl_range` is
    more than `threshold`x (or less than 1/`threshold`x) the group's median -
    e.g. two batches of the same "site reflectance" file calibrated/processed
    differently. Harmless for ratio-only algorithms (scale cancels), but
    band-ratio-plus-product equations like the MODIS one aren't fully
    scale-invariant, so a mixed-scale input site set will show up as spurious
    extra scatter/bias that isn't real signal.
    """
    wl_cols = [c for c in df.columns if isinstance(c, (int, float)) and wl_range[0] <= c <= wl_range[1]]
    if not wl_cols:
        return
    row_means = df[wl_cols].mean(axis=1).to_numpy(dtype=float)
    median = np.median(row_means)
    if median == 0:
        return
    ratios = row_means / median
    names = df["wavelength"] if "wavelength" in df.columns else df.index
    flagged = [(name, ratio) for name, ratio in zip(names, ratios) if ratio > threshold or ratio < 1 / threshold]
    if flagged:
        print(
            f"WARNING: {len(flagged)}/{len(df)} site(s) have mean reflectance "
            f">{threshold:g}x (or <1/{threshold:g}x) the group median over "
            f"{wl_range[0]}-{wl_range[1]}nm - likely a different processing "
            f"batch/calibration mixed into this file. Not fully harmless for "
            f"non-ratio-only algorithms (e.g. MODIS): {flagged}"
        )


def load_dataset_via_srf(path: str, bands_nm: "tuple[int, ...]", compute_x, srf_path: str, sheet_name: Optional[str] = None) -> "tuple[np.ndarray, np.ndarray]":
    df = _read_table(path, sheet_name)
    flag_scale_outliers(df)
    modis_df = convolve_to_modis_bands(df, srf_path)
    return _x_and_y_from_df(modis_df, bands_nm, compute_x)


def get_model_errors(actual: np.ndarray, predicted: np.ndarray) -> "tuple[float, float, float, float, float, float]":
    """RMSE/MAE in mg/L and as a %-of-mean, plus Bias and rRMSE%.

    RMSE(%)/MAE(%) normalize by the mean of `actual` (how big the error is
    relative to the typical SPM value); rRMSE(%) instead averages each point's
    *own* relative error before taking the root - a different, commonly
    reported statistic that penalizes errors on low-SPM samples more heavily.
    """
    actual = np.asarray(actual, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    diff = predicted - actual

    rmse = float(np.sqrt(np.mean(diff ** 2)))
    mae = float(np.mean(np.abs(diff)))
    bias = float(np.mean(diff))

    mean_actual = float(np.mean(actual))
    rmse_pct = rmse / mean_actual * 100
    mae_pct = mae / mean_actual * 100
    rrmse_pct = float(np.sqrt(np.mean((diff / actual) ** 2)) * 100)

    return rmse, rmse_pct, mae, mae_pct, bias, rrmse_pct


def full_model_stats(x: np.ndarray, y: np.ndarray) -> "tuple[float, float, int, float]":
    """R^2, F-statistic, residual degrees of freedom, and significance (p-value)
    for the simple OLS fit of y ~ x on the *whole* training set - the "Table 3"
    style summary stats, as opposed to the k-fold numbers below.

    With one predictor, the two-sided p-value scipy.stats.linregress reports for
    the slope is identical to the overall F-test's p-value, so it doubles as
    "Significance" here without a separate F-distribution CDF computation.
    """
    n = len(x)
    slope, intercept, r_value, p_value, _std_err = scipy_stats.linregress(x, y)
    r2 = r_value ** 2
    pred = slope * x + intercept
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    df_resid = n - 2
    f_stat = ((ss_tot - ss_res) / 1) / (ss_res / df_resid)
    return float(r2), float(f_stat), int(df_resid), float(p_value)


def cross_validate_errors(x: np.ndarray, y: np.ndarray, n_splits: int = 10, random_state: int = 42) -> Dict[str, np.ndarray]:
    """Refit on each of `n_splits` folds of the training set and score on the
    held-out fold, mirroring this project's original dissertation methodology
    (K_hyperspectral_total_rrs-Error_Result.csv: 10 folds of the 50-sample
    training set) - this is where Table 3's RMSE/MAE come from, not panel (b)'s
    external validation set.
    """
    n_splits = min(n_splits, len(x))
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    per_fold = {"rmse": [], "rmse_pct": [], "mae": [], "mae_pct": []}
    for train_idx, test_idx in kf.split(x):
        fold_model = LinearRegression()
        fold_model.fit(x[train_idx].reshape(-1, 1), y[train_idx])
        pred = fold_model.predict(x[test_idx].reshape(-1, 1))
        rmse, rmse_pct, mae, mae_pct, _bias, _rrmse_pct = get_model_errors(y[test_idx], pred)
        per_fold["rmse"].append(rmse)
        per_fold["rmse_pct"].append(rmse_pct)
        per_fold["mae"].append(mae)
        per_fold["mae_pct"].append(mae_pct)
    return {key: np.array(values) for key, values in per_fold.items()}


def write_statistics_csv(path: str, r2: float, f_stat: float, df_resid: int, p_value: float, cv: Dict[str, np.ndarray]) -> None:
    # "Significance" here is displayed like a paper's table (e.g. "<0.001"), not
    # the raw float - the exact p-value from a near-perfect single-predictor fit
    # (e.g. 7e-15) isn't meaningfully different from any other p < 0.001 for
    # reporting purposes.
    significance = "<0.001" if p_value < 0.001 else f"{p_value:.3f}"
    rows = [
        ("R-squared", round(r2, 2), "", ""),
        ("F-statistic", round(f_stat, 1), "", ""),
        ("Degree of Freedom", df_resid, "", ""),
        ("Significance", significance, "", ""),
        ("Root Mean Squared Error (RMSE) (mg/L)", round(cv["rmse"].mean(), 2), round(cv["rmse"].min(), 2), round(cv["rmse"].max(), 2)),
        ("Root Mean Squared Error (RMSE) (percent)", round(cv["rmse_pct"].mean(), 1), round(cv["rmse_pct"].min(), 1), round(cv["rmse_pct"].max(), 1)),
        ("Mean Absolute Error (MAE) (mg/L)", round(cv["mae"].mean(), 2), round(cv["mae"].min(), 2), round(cv["mae"].max(), 2)),
        ("Mean Absolute Error (MAE) (percent)", round(cv["mae_pct"].mean(), 1), round(cv["mae_pct"].min(), 1), round(cv["mae_pct"].max(), 1)),
    ]
    pd.DataFrame(rows, columns=["Statistical Results", "Value", "Min", "Max"]).to_csv(path, index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--algorithm", choices=sorted(ALGORITHMS), default="hyperspectral", help="Which published SPM algorithm to validate.")
    parser.add_argument("--training-path", default=None, help="CSV/Excel with the algorithm's Rrs bands + SPM columns, used for panel (a)'s diagnostic scatter.")
    parser.add_argument("--testing-path", default=None, help="Excel/CSV with the held-out in-situ validation data for panel (b).")
    parser.add_argument("--testing-sheet", default=None, help="Sheet name in --testing-path (Excel only).")
    parser.add_argument("--testing-srf", default=None, help="Spectral response function CSV; if set (or the algorithm has a default one), --testing-path is treated as per-nm hyperspectral Rrs and band-convolved down to this algorithm's bands first.")
    parser.add_argument("--output", default=None, help="Output figure path; defaults to data/output/<algorithm>_algorithm_output.svg")
    args = parser.parse_args()

    algo = ALGORITHMS[args.algorithm]
    training_path = args.training_path or algo["default_training_path"]
    testing_path = args.testing_path or algo["default_testing_path"]
    testing_sheet = args.testing_sheet or algo["default_testing_sheet"]
    testing_srf_path = args.testing_srf or algo.get("testing_srf_path")

    if training_path is None or testing_path is None:
        sys.exit(
            f"No bundled training/testing data for --algorithm {args.algorithm} yet. "
            f"Pass --training-path and --testing-path explicitly."
        )

    output_dir = os.path.join(REPO_ROOT, "data", "output")
    os.makedirs(output_dir, exist_ok=True)
    output_path = args.output or os.path.join(output_dir, f"{args.algorithm}_algorithm_output.svg")
    stats_csv_path = os.path.join(output_dir, f"{args.algorithm}_algorithm_statistics.csv")

    x_train, y_train = load_dataset(training_path, algo["bands_nm"], algo["compute_x"], x_column=algo.get("default_training_x_column"))
    local_model = LinearRegression()
    local_model.fit(x_train.reshape(-1, 1), y_train)
    r2_train = local_model.score(x_train.reshape(-1, 1), y_train)

    if testing_srf_path:
        x_test, y_test = load_dataset_via_srf(testing_path, algo["bands_nm"], algo["compute_x"], testing_srf_path, sheet_name=testing_sheet)
        testing_bands_csv_path = os.path.join(output_dir, f"{args.algorithm}_algorithm_testing_bands.csv")
        raw_testing_df = _read_table(testing_path, testing_sheet)
        build_modis_band_table(raw_testing_df, testing_srf_path).to_csv(testing_bands_csv_path, index=False)
        print(f"Saved testing SPM/MODIS-band table to {testing_bands_csv_path}")
    else:
        x_test, y_test = load_dataset(testing_path, algo["bands_nm"], algo["compute_x"], sheet_name=testing_sheet)
    fixed_model = FixedLinearModel(algo["slope"], algo["intercept"])
    predicted = fixed_model.predict(x_test)
    actual = y_test

    rmse, rmse_pct, mae, mae_pct, bias, rrmse_pct = get_model_errors(actual, predicted)
    validation = OrderedDict([
        ("RMSE(mg/L)", rmse),
        ("RMSE(%)", rmse_pct),
        ("MAE(mg/L)", mae),
        ("MAE(%)", mae_pct),
        ("Bias(mg/L)", bias),
        ("rRMSE(%)", rrmse_pct),
    ])

    print(f"Algorithm: {args.algorithm}  (published: SPM = {algo['slope']}*X + {algo['intercept']})")
    print(f"Training ({training_path}): n={len(x_train)}")
    print(f"  local best fit: y = {local_model.coef_[0]:.2f}x + {local_model.intercept_:.2f}   R^2={r2_train:.3f}")
    print(f"Testing  ({testing_path}{' [' + testing_sheet + ']' if testing_sheet else ''}): n={len(x_test)}")
    for key, value in validation.items():
        print(f"  {key}: {value:.2f}")

    r2_full, f_stat, df_resid, p_value = full_model_stats(x_train, y_train)
    cv = cross_validate_errors(x_train, y_train)
    write_statistics_csv(stats_csv_path, r2_full, f_stat, df_resid, p_value, cv)
    print(f"R-squared: {r2_full:.2f}  F-statistic: {f_stat:.1f}  Degree of Freedom: {df_resid}  Significance: {p_value:.3f}")
    print(f"  RMSE(mg/L) mean/min/max: {cv['rmse'].mean():.2f} / {cv['rmse'].min():.2f} / {cv['rmse'].max():.2f}")
    print(f"  MAE(mg/L)  mean/min/max: {cv['mae'].mean():.2f} / {cv['mae'].min():.2f} / {cv['mae'].max():.2f}")
    print(f"Saved statistics table to {stats_csv_path}")

    plot.create_scatter_plot(
        list(actual), list(predicted), validation, local_model,
        list(y_train), list(x_train), r2_train, algo["equation_latex"],
        size="paper", output_path=output_path,
    )
    print(f"Saved figure to {output_path}")


if __name__ == "__main__":
    main()
