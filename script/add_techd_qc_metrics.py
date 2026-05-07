#!/usr/bin/env python3
import argparse
import os
import sys
import time
from typing import List, Dict

import math

import numpy as np
import pandas as pd
import yaml

# -------------------------
# Config / constants
# -------------------------
EM_BOUNDARY = 1.05
ML_BOUNDARY = 1.10
RT_PHASE_MAP = ["none", "early", "late", "peak"]

NEW_YAML_COLUMNS: List[Dict[str, str]] = [
    {"dtype": "float", "name": "reads_meanlog"},
    {"dtype": "float", "name": "reads_sdlog"},
    {"dtype": "float", "name": "reads_gc_cor"},
    {"dtype": "float", "name": "corrected_mean"}, # likely important
    {"dtype": "float", "name": "corrected_sd"}, # likely important
    {"dtype": "float", "name": "corrected_bottom"}, # important
    {"dtype": "float", "name": "corrected_diploid_perc"}, # important
    {"dtype": "float", "name": "rt_early"},
    {"dtype": "float", "name": "rt_mid"},
    {"dtype": "float", "name": "rt_late"},
    {"dtype": "float", "name": "rt_el"},
    {"dtype": "float", "name": "rt_em"},
    {"dtype": "float", "name": "rt_ml"},
    {"dtype": "float", "name": "chr_sd"}, # important
    {"dtype": "float", "name": "rt_em_split"},
    {"dtype": "float", "name": "rt_ml_split"},
    {"dtype": "float", "name": "rt_phase_code"},
    {"dtype": "str",   "name": "rt_phase_name"}, # important
]


def parse_args():
    p = argparse.ArgumentParser(
        description="Add technical depth QC metrics to metrics CSV.gz and append schema to YAML."
    )
    p.add_argument("--reads-csv", required=True,
                   help="Input reads CSV.gz (e.g., <cell_id>_reads.csv.gz)")
    p.add_argument("--metrics-csv", required=True,
                   help="Input/output metrics CSV.gz (e.g., <cell_id>_metrics.csv.gz)")
    p.add_argument("--metrics-yaml", required=True,
                   help="Input/output YAML schema (e.g., <cell_id>_metrics.csv.gz.yaml)")
    p.add_argument("--htert-rt", default="htert_rt.csv",
                   help="Path to htert_rt.csv replication-timing reference (default: htert_rt.csv)")
    p.add_argument("--em-boundary", type=float, default=EM_BOUNDARY,
                   help=f"rt_em threshold for S-phase early/peak split (default: {EM_BOUNDARY})")
    p.add_argument("--ml-boundary", type=float, default=ML_BOUNDARY,
                   help=f"rt_ml threshold for S-phase late/peak split (default: {ML_BOUNDARY})")
    p.add_argument("--min-median-reads", type=float, default=5.0,
                   help="Minimum median read depth per bin to process a cell (default: 5)")
    return p.parse_args()


def load_yaml(path: str) -> dict:
    with open(path, "r") as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict):
        raise ValueError(f"{path} must be a YAML mapping at top level")
    data.setdefault("columns", [])
    return data


def save_yaml(doc: dict, path: str):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        yaml.safe_dump(doc, fh, sort_keys=False, default_flow_style=False, allow_unicode=True)
    os.replace(tmp, path)


def ensure_yaml_columns(doc: dict, new_cols: List[Dict[str, str]]):
    existing = {c.get("name") for c in doc.get("columns", []) if isinstance(c, dict)}
    for c in new_cols:
        if c["name"] not in existing:
            doc["columns"].append({"dtype": c["dtype"], "name": c["name"]})
    doc.setdefault("header", True)
    doc.setdefault("sep", ",")


def _nan_result() -> dict:
    result = {}
    for col in NEW_YAML_COLUMNS:
        result[col["name"]] = np.nan
    return result


def process_cell(cell_id: str, cell0: pd.DataFrame, htert: pd.DataFrame,
                 min_median_reads: float) -> dict:
    """Compute per-cell techd metrics. Returns a dict of new column values (NaN for skipped cells)."""
    result = _nan_result()

    if cell0.empty or cell0["reads"].median() < min_median_reads:
        print(f"  {cell_id}: too few reads, skipping", file=sys.stderr)
        return result

    # Log-normal fit on raw reads
    log_reads = np.log(cell0["reads"].astype(float))
    meanlog = log_reads.mean()
    sdlog = log_reads.std(ddof=1)
    cell0 = cell0.copy()
    cell0["z"] = (log_reads - meanlog) / sdlog
    cell0["outlier"] = cell0["z"].abs() > 3

    # Pearson correlation between gc and reads (matches R's cor.test default)
    r_val = float(np.corrcoef(cell0["gc"].astype(float), cell0["reads"].astype(float))[0, 1])

    result["reads_meanlog"] = round(float(meanlog), 3)
    result["reads_sdlog"] = round(float(sdlog), 3)
    result["reads_gc_cor"] = round(float(r_val), 3)

    # S-phase metrics — join on coordinates only to avoid float precision issues with gc/map
    if htert is not None and not htert.empty:
        coord_cols = [c for c in ("chr", "start", "end") if c in cell0.columns and c in htert.columns]
        cell_s = cell0.merge(htert[coord_cols + ["htertrt"]], on=coord_cols, how="inner")

        if cell_s.empty:
            print(f"  {cell_id}: hTERT reference unable to merge, likely not hg19 reference",
                  file=sys.stderr)
        else:
            rt_stats = (
                cell_s[cell_s["htertrt"].notna()]
                .groupby("htertrt", observed=True)["reads"]
                .mean()
                .rename_axis("htertrt")
                .reset_index(name="mreads")
            )
            rt = dict(zip(rt_stats["htertrt"], rt_stats["mreads"]))

            for phase in ("early", "mid", "late"):
                if phase in rt:
                    result[f"rt_{phase}"] = float(rt[phase])

            early, mid, late = rt.get("early"), rt.get("mid"), rt.get("late")
            if early is not None and mid and mid != 0:
                result["rt_em"] = round(early / mid, 3)
            if early is not None and late and late != 0:
                result["rt_el"] = round(early / late, 3)
            if mid is not None and late and late != 0:
                result["rt_ml"] = round(mid / late, 3)

    # Normalized read variance on non-outlier bins
    cell = cell0[~cell0["outlier"]]

    if "cor_gc" not in cell.columns or cell.empty:
        return result

    joint_mean = float(cell["cor_gc"].mean())
    joint_sd = float(cell["cor_gc"].std(ddof=1))

    if joint_sd > 0:
        _norm_cdf = lambda x: 0.5 * (1.0 + math.erf((x - joint_mean) / (joint_sd * math.sqrt(2))))
        range_perc = _norm_cdf(1.5) - _norm_cdf(0.5)
    else:
        range_perc = np.nan

    result["corrected_mean"] = round(joint_mean, 3)
    result["corrected_sd"] = round(joint_sd, 3)
    result["corrected_bottom"] = round(joint_mean - 2 * joint_sd, 3)
    result["corrected_diploid_perc"] = round(range_perc, 3)

    # Chromosomal aberration: SD of per-chromosome mean cor_gc
    chr_means = cell.groupby("chr")["cor_gc"].mean()
    result["chr_sd"] = round(float(chr_means.std(ddof=1)), 3)

    return result


def categorize_sphase(metrics: pd.DataFrame, em_boundary: float, ml_boundary: float) -> pd.DataFrame:
    """Add rt_em_split, rt_ml_split, rt_phase_code, rt_phase_name; propagates NaN."""
    em = metrics["rt_em"]
    ml = metrics["rt_ml"]

    metrics["rt_em_split"] = em.where(em.isna(), em >= em_boundary)
    metrics["rt_ml_split"] = ml.where(ml.isna(), ml >= ml_boundary)

    phase_code = (
        metrics["rt_em_split"].astype(float) * 1 +
        metrics["rt_ml_split"].astype(float) * 2
    )
    metrics["rt_phase_code"] = phase_code.where(phase_code.notna(), other=np.nan)
    metrics["rt_phase_name"] = phase_code.apply(
        lambda x: RT_PHASE_MAP[int(x)] if pd.notna(x) else None
    )
    return metrics


def main():
    args = parse_args()
    start = time.time()

    htert = None
    if not os.path.exists(args.htert_rt):
        print(f"Warning: {args.htert_rt} not found — S-phase metrics will be skipped",
              file=sys.stderr)
    else:
        htert = pd.read_csv(args.htert_rt)

    print(f"{os.path.basename(args.metrics_csv)}: starting process", file=sys.stderr)

    metrics = pd.read_csv(args.metrics_csv)
    reads = pd.read_csv(args.reads_csv)

    # initialize new columns with NaN
    for col in NEW_YAML_COLUMNS:
        metrics[col["name"]] = np.nan

    # pre-filter and group reads by cell to avoid repeated full-table scans
    reads_filtered = reads[(reads["reads"] != 0) & (reads["gc"] != -1)]
    reads_by_cell = {cid: grp for cid, grp in reads_filtered.groupby("cell_id")}

    for i, row in metrics.iterrows():
        cell_id = row["cell_id"]
        cell0 = reads_by_cell.get(cell_id, pd.DataFrame())
        result = process_cell(cell_id, cell0, htert, args.min_median_reads)
        for key, val in result.items():
            metrics.at[i, key] = val

    metrics = categorize_sphase(metrics, args.em_boundary, args.ml_boundary)

    elapsed = round(time.time() - start, 2)
    print(f"{os.path.basename(args.metrics_csv)}: total time: {elapsed} seconds", file=sys.stderr)

    # write CSV back in place (compressed)
    tmp_csv = args.metrics_csv + ".tmp"
    metrics.to_csv(tmp_csv, index=False, compression="gzip")
    os.replace(tmp_csv, args.metrics_csv)

    # update YAML schema
    schema = load_yaml(args.metrics_yaml)
    ensure_yaml_columns(schema, NEW_YAML_COLUMNS)
    save_yaml(schema, args.metrics_yaml)

    print(f"Updated CSV:  {args.metrics_csv}")
    print(f"Updated YAML: {args.metrics_yaml}")


if __name__ == "__main__":
    # deps: pandas, pyyaml, numpy
    main()
