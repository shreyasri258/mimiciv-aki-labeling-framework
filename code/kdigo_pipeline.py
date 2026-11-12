"""
Shared, parameterized core functions for KDIGO AKI labeling.
This is used by both AKI labeling & outcome prevalance notebook and the sensitivity-analysis notebook.
"""

import numpy as np
import pandas as pd


# Measurement-gap computation
def compute_measurement_gaps(outputevents_df):
    """Compute minutes elapsed since the previous measurement, per stay."""
    temp = outputevents_df.sort_values(["stay_id", "charttime"]).copy()
    temp["tm_since_last"] = (
        temp.groupby("stay_id")["charttime"].diff().dt.total_seconds().fillna(0) / 60
    )
    return temp


# Urine-output criterion (elapsed-time method)
def kdigo_urine_elapsed_time(stay_events, weight_kg, window_start, window_end, max_gap_hours=None):
    """
    Evaluate urine-output evaluability for one stay using the elapsed-time method.
    """
    if pd.isna(weight_kg) or weight_kg <= 0:
        return False, False

    stay_events = stay_events.sort_values("charttime").reset_index(drop=True)
    any_evaluable = False
    any_positive = False

    for i, row in stay_events.iterrows():
        # anchor must fall inside the target window to prevent label leakage
        if not (window_start <= row["hours_from_icu"] < window_end):
            continue

        window = stay_events[
            (stay_events["charttime"] <= row["charttime"]) &
            (stay_events["charttime"] >= row["charttime"] - pd.Timedelta(hours=6))
        ]

        if max_gap_hours is not None:
            max_gap_in_window = window["tm_since_last"].max() / 60 if len(window) else 0
            if max_gap_in_window > max_gap_hours:
                continue  # window contains a gap too large to trust; skip this anchor

        elapsed_hours = window["tm_since_last"].sum() / 60
        if elapsed_hours >= 6:
            any_evaluable = True
            rate = window["urine_ml"].sum() / weight_kg / elapsed_hours
            if rate < 0.5:
                any_positive = True

    return any_evaluable, any_positive


def generate_urine_labels_elapsed_time(full_stay_events_df, weight_df, window_start, window_end, max_gap_hours=None):
    """Compute per-stay urine-output evaluability and AKI-positivity for a full cohort."""
    weight_map = weight_df.set_index("stay_id")["weight_kg"]
    results = []
    for stay_id, group in full_stay_events_df.groupby("stay_id"):
        wt = weight_map.get(stay_id, np.nan)
        evaluable, positive = kdigo_urine_elapsed_time(
            group, wt, window_start, window_end, max_gap_hours=max_gap_hours
        )
        results.append({"stay_id": stay_id, "urine_evaluable_et": evaluable, "urine_labels_et": positive})
    return pd.DataFrame(results)



# Serum creatinine criteria
def kdigo_scr_15x(creatinine_df, baseline_df):
    """KDIGO Criterion: Serum creatinine >= 1.5 x baseline."""
    temp = creatinine_df.copy()
    if "baseline_scr" not in temp.columns:
        temp = temp.merge(baseline_df[["stay_id", "baseline_scr"]], on="stay_id", how="left")
    temp["aki_1_5x"] = (temp["valuenum"] >= 1.5 * temp["baseline_scr"])
    return temp.groupby("stay_id")["aki_1_5x"].max().reset_index()


def kdigo_scr_03(creatinine_df, outcome_start=24, outcome_end=48, lookback_hrs=48):
    """KDIGO Criterion: Serum creatinine rise >= 0.3 mg/dL within lookback_hrs,
    with the later qualifying measurement inside [outcome_start, outcome_end)."""
    temp = creatinine_df[["stay_id", "charttime", "valuenum", "hours_from_icu"]].copy()
    temp = temp.dropna(subset=["valuenum", "charttime"]).sort_values(["stay_id", "charttime"])

    later = temp[(temp["hours_from_icu"] >= outcome_start) & (temp["hours_from_icu"] < outcome_end)].copy()
    earlier = temp[temp["hours_from_icu"] >= (outcome_start - lookback_hrs)].copy()

    pairs = earlier.merge(later, on="stay_id", suffixes=("_earlier", "_later"))
    pairs["hours_diff"] = (pairs["charttime_later"] - pairs["charttime_earlier"]).dt.total_seconds() / 3600
    pairs = pairs[(pairs["hours_diff"] > 0) & (pairs["hours_diff"] <= lookback_hrs)].copy()
    pairs["scr_change"] = (pairs["valuenum_later"] - pairs["valuenum_earlier"]).round(2)
    pairs["aki_0_3"] = pairs["scr_change"] >= 0.3

    return pairs.groupby("stay_id")["aki_0_3"].max().reset_index()


def kdigo_scr_03_early(creatinine_df):
    """Early-AKI assessment using only 0-24h measurements (prevents outcome-window leakage)."""
    temp = creatinine_df[["stay_id", "charttime", "valuenum"]].copy()
    temp = temp.dropna(subset=["valuenum"]).sort_values(["stay_id", "charttime"])
    pairs = temp.merge(temp, on="stay_id", suffixes=("_earlier", "_later"))
    pairs["hours_diff"] = (pairs["charttime_later"] - pairs["charttime_earlier"]).dt.total_seconds() / 3600
    pairs = pairs[(pairs["hours_diff"] > 0) & (pairs["hours_diff"] <= 24)].copy()
    pairs["scr_change"] = (pairs["valuenum_later"] - pairs["valuenum_earlier"]).round(2)
    pairs["aki_0_3"] = pairs["scr_change"] >= 0.3
    return pairs.groupby("stay_id")["aki_0_3"].max().reset_index()


def combine_scr_labels(cohort_df, aki15_df, aki03_df):
    """Combine the two KDIGO serum-creatinine sub-criteria into one aki_scr flag."""
    labels = cohort_df[["stay_id"]].merge(aki15_df, on="stay_id", how="left")
    labels = labels.merge(aki03_df, on="stay_id", how="left")
    labels["aki_1_5x"] = labels["aki_1_5x"].fillna(False)
    labels["aki_0_3"] = labels["aki_0_3"].fillna(False)
    labels["aki_scr"] = labels["aki_1_5x"] | labels["aki_0_3"]
    return labels


def get_time_window(df, start_hour, end_hour):
    """Return rows within a half-open ICU time window [start_hour, end_hour)."""
    return df[(df["hours_from_icu"] >= start_hour) & (df["hours_from_icu"] < end_hour)].copy()



# Baseline creatinine 
def compute_baseline_scr(creatinine_df, cohort_df, method="fixed"):
    """
    Compute per-stay baseline creatinine.
    method="fixed" (default behavior):
    method="rolling" (variant):
       
    """
    if method == "fixed":
        baseline_window = creatinine_df[
            (creatinine_df["hours_from_icu"] >= -168) & (creatinine_df["hours_from_icu"] < 0)
        ]
        baseline_scr = (
            baseline_window.groupby("stay_id")["valuenum"].min()
            .reset_index().rename(columns={"valuenum": "baseline_scr"})
        )
        missing_baseline = cohort_df[["stay_id"]].merge(baseline_scr, on="stay_id", how="left")

        first_scr = (
            creatinine_df[creatinine_df["hours_from_icu"] >= 0]
            .sort_values("charttime").groupby("stay_id").first()[["valuenum"]]
            .rename(columns={"valuenum": "first_scr"}).reset_index()
        )
        baseline_scr = missing_baseline.merge(first_scr, on="stay_id", how="left")
        baseline_scr["pre_icu_baseline"] = baseline_scr["baseline_scr"]
        baseline_scr["baseline_scr"] = baseline_scr["baseline_scr"].fillna(baseline_scr["first_scr"])
        return baseline_scr

    elif method == "rolling":
        compute_rolling_creatinine_baseline(creatinine_df)

    else:
        raise ValueError(f"Unknown baseline method: {method!r}")



# Urine-output evaluability handling (variant)
def get_urine_evaluability_labels(aki_labels, mode="exclude"):
    """
    mode="exclude": stays with no baseline AND no urine evaluability
        are dropped upstream in cohort construction.
    mode="negative": urine-unevaluable stays kept, urine_labels forced False.
    mode="strict_exclude_unevaluable": ALL urine-unevaluable stays removed from the cohort.
    """
    if mode == "exclude":
        return aki_labels
    elif mode == "negative":
        aki_labels = aki_labels.copy()
        aki_labels.loc[~aki_labels["urine_evaluable"], "urine_labels"] = False
        aki_labels["AKI"] = aki_labels["aki_scr"] | aki_labels["urine_labels"]
        return aki_labels
    elif mode == "strict_exclude_unevaluable":
        return aki_labels[aki_labels["urine_evaluable"]].copy()
    else:
        raise ValueError(f"Unknown evaluability mode: {mode!r}")


# Weight determination (variant)

def compute_weight_table(chartevents_df, method="first_24h"):
    """
    Compute per-stay usable weight.

    method="first_24h" (default):
        First kg-native (or lb-converted, kg-priority) weight charted within
        the first 24h of the ICU stay. Held fixed for the whole stay.

    method="full_stay_last_known":
        searches the ENTIRE stay (not just first 24h) for the first available weight. Tests whether the
        24h time restriction itself (not data source) costs usable coverage.
    """
    WEIGHT_ITEMID_KG = 226512
    WEIGHT_ITEMID_KG_DAILY = 224639
    WEIGHT_ITEMID_LB = 226531

    w_raw = chartevents_df[
        chartevents_df["itemid"].isin([WEIGHT_ITEMID_KG, WEIGHT_ITEMID_KG_DAILY, WEIGHT_ITEMID_LB])
    ].copy()
    w_raw["hours_from_icu"] = (w_raw["charttime"] - w_raw["intime"]).dt.total_seconds() / 3600

    if method == "first_24h":
        w_raw = w_raw[(w_raw["hours_from_icu"] >= 0) & (w_raw["hours_from_icu"] < 24)].copy()
    elif method == "full_stay_last_known":
        w_raw = w_raw[w_raw["hours_from_icu"] >= 0].copy()
    else:
        raise ValueError(f"Unknown weight method: {method!r}")

    w_raw["weight_kg"] = np.where(
        w_raw["itemid"] == WEIGHT_ITEMID_LB,
        pd.to_numeric(w_raw["valuenum"], errors="coerce") * 0.45359237,
        pd.to_numeric(w_raw["valuenum"], errors="coerce")
    )
    w_raw = w_raw.dropna(subset=["weight_kg"])
    w_raw = w_raw[w_raw["weight_kg"].between(20, 300)].copy()
    w_raw["priority"] = np.where(
        w_raw["itemid"].isin([WEIGHT_ITEMID_KG, WEIGHT_ITEMID_KG_DAILY]), 0, 1
    )
    w_raw = w_raw.sort_values(["stay_id", "priority", "charttime"])
    weight = w_raw.groupby("stay_id", as_index=False).first()[["stay_id", "weight_kg"]]
    return weight


def compute_weight_from_durations(weight_durations_df, cohort_df):
    wd = weight_durations_df.merge(cohort_df[["stay_id"]], on="stay_id", how="inner")
    wd = wd.sort_values(["stay_id", "starttime"])
    weight = wd.groupby("stay_id", as_index=False).first()[["stay_id", "weight"]]
    weight = weight.rename(columns={"weight": "weight_kg"})
    weight = weight[weight["weight_kg"].between(20, 300)].copy()
    return weight


def compute_rolling_creatinine_baseline(creatinine_df):
    temp = creatinine_df[["stay_id", "charttime", "valuenum"]].copy()
    temp = temp.dropna(subset=["valuenum"]).sort_values(["stay_id", "charttime"]).reset_index(drop=True)

    merged = temp.merge(temp, on="stay_id", suffixes=("", "_prior"))
    merged = merged[merged["charttime_prior"] < merged["charttime"]]  # CHANGED: strictly before
    merged["hours_diff"] = (merged["charttime"] - merged["charttime_prior"]).dt.total_seconds() / 3600

    past_48h = merged[merged["hours_diff"] <= 48]
    creat_low_48 = (
        past_48h.groupby(["stay_id", "charttime"])["valuenum_prior"].min()
        .reset_index().rename(columns={"valuenum_prior": "creat_low_past_48hr"}))

    past_7d = merged[merged["hours_diff"] <= 168]
    creat_low_7d = (
        past_7d.groupby(["stay_id", "charttime"])["valuenum_prior"].min()
        .reset_index().rename(columns={"valuenum_prior": "creat_low_past_7day"})
    )

    result = temp.merge(creat_low_48, on=["stay_id", "charttime"], how="left")
    result = result.merge(creat_low_7d, on=["stay_id", "charttime"], how="left")
    return result
def kdigo_scr_15x_rolling(creatinine_rolling_df, outcome_start=None, outcome_end=None, hours_from_icu_df=None):
    """
    Rolling-baseline version of the 1.5x criterion: compares each measurement against its own trailing 7-day minimum, rather than a single fixed baseline.
    If outcome_start/outcome_end given, restricts to measurements with hours_from_icu in that window."""
    temp = creatinine_rolling_df.copy()
    if hours_from_icu_df is not None:
        temp = temp.merge(hours_from_icu_df[["stay_id", "charttime", "hours_from_icu"]], on=["stay_id", "charttime"], how="left")
    if outcome_start is not None and outcome_end is not None:
        temp = temp[(temp["hours_from_icu"] >= outcome_start) & (temp["hours_from_icu"] < outcome_end)]
    temp["aki_1_5x"] = temp["valuenum"] >= 1.5 * temp["creat_low_past_7day"]
    return temp.groupby("stay_id")["aki_1_5x"].max().reset_index()


def kdigo_scr_03_rolling(creatinine_rolling_df, outcome_start=24, outcome_end=48, hours_from_icu_df=None):
    """
    Rolling-baseline version of the 0.3 mg/dL criterion: compares each measurement in the outcome window against ITS own trailing 48h minimum.
    """
    temp = creatinine_rolling_df.copy()
    if hours_from_icu_df is not None:
        temp = temp.merge(hours_from_icu_df[["stay_id", "charttime", "hours_from_icu"]], on=["stay_id", "charttime"], how="left")
    temp = temp[(temp["hours_from_icu"] >= outcome_start) & (temp["hours_from_icu"] < outcome_end)]
    temp["aki_0_3"] = (temp["valuenum"] - temp["creat_low_past_48hr"]).round(2) >= 0.3
    return temp.groupby("stay_id")["aki_0_3"].max().reset_index()

