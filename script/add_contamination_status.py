#!/usr/bin/env python3
import argparse
import os
import sys
from typing import List, Dict

import pandas as pd
import yaml

# -------------------------
# Config / constants
# -------------------------
NEW_YAML_COLUMNS: List[Dict[str, str]] = [
    {"dtype": "float", "name": "fastqscreen_eukaryotes"},
    {"dtype": "float", "name": "fastqscreen_prokaryotes"},
    {"dtype": "float", "name": "fastqscreen_sum"},
    {"dtype": "float", "name": "contam_eukar_perc"},
    {"dtype": "float", "name": "contam_human_perc"},
    {"dtype": "float", "name": "contam_mouse_perc"},
    {"dtype": "float", "name": "contam_salmo_perc"},
    {"dtype": "float", "name": "contam_nohit_perc"},
    {"dtype": "float", "name": "contam_pseud_perc"},
    {"dtype": "str",   "name": "contam_species"},
    {"dtype": "str",   "name": "contam_species_purity"},
]

REQUIRED_COLS = [
    "fastqscreen_human",
    "fastqscreen_mouse",
    "fastqscreen_salmon",
    "fastqscreen_nohit",
    "fastqscreen_pseudomonas",
]


def parse_args():
    p = argparse.ArgumentParser(
        description="Add contamination-derived columns to metrics CSV.gz and append schema to YAML."
    )
    p.add_argument("--metrics-csv", required=True,
                   help="Input metrics CSV.gz (e.g., <cell_id>_metrics.csv.gz)")
    p.add_argument("--metrics-yaml", required=True,
                   help="Input YAML schema (e.g., <cell_id>_metrics.csv.gz.yaml)")
    # kept for backward compat with your current call; if given, it overwrites --metrics-csv
    p.add_argument("--input", help="(legacy) same as --metrics-csv")
    p.add_argument("--output", help="(legacy) ignored; we overwrite in place")
    p.add_argument("--org-threshold", type=float, default=0.60,
                   help="Threshold for HUMAN/MOUSE call (default: 0.60)")
    p.add_argument("--euk-fraction-clean-min", type=float, default=0.65,
                   help="contam_eukar_perc < this => CONTAMINATED (default: 0.65)")
    p.add_argument("--euk-fraction-negative-min", type=float, default=0.30,
                   help="contam_eukar_perc < this => NEGATIVE (default: 0.30)")
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
    # preserve header/sep defaults if missing
    doc.setdefault("header", True)
    # the YAML you showed has sep: ',' (a string with single quotes in YAML),
    # dumping as ',' is fine; readers treat it the same.
    doc.setdefault("sep", ",")


def compute_derived(df: pd.DataFrame, org_threshold: float,
                    clean_min: float, negative_min: float) -> pd.DataFrame:
    # validate required columns
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError("Wrong version of metrics; missing: " + ", ".join(missing))

    # sums
    df["fastqscreen_eukaryotes"] = (
        df["fastqscreen_human"].astype(float) +
        df["fastqscreen_mouse"].astype(float) +
        df["fastqscreen_salmon"].astype(float)
    )
    df["fastqscreen_prokaryotes"] = (
        df["fastqscreen_nohit"].astype(float) +
        df["fastqscreen_pseudomonas"].astype(float)
    )
    df["fastqscreen_sum"] = df["fastqscreen_eukaryotes"] + df["fastqscreen_prokaryotes"]

    # denominators with NA to avoid div-by-zero
    denom_sum = df["fastqscreen_sum"].replace(0, pd.NA)
    denom_euk = df["fastqscreen_eukaryotes"].replace(0, pd.NA)
    denom_pro = df["fastqscreen_prokaryotes"].replace(0, pd.NA)

    # fractions
    df["contam_eukar_perc"] = df["fastqscreen_eukaryotes"] / denom_sum
    df["contam_human_perc"] = df["fastqscreen_human"] / denom_euk
    df["contam_mouse_perc"] = df["fastqscreen_mouse"] / denom_euk
    df["contam_salmo_perc"] = df["fastqscreen_salmon"] / denom_euk
    df["contam_nohit_perc"] = df["fastqscreen_nohit"] / denom_pro
    df["contam_pseud_perc"] = df["fastqscreen_pseudomonas"] / denom_pro

    # classification
    df["contam_species"] = "MIXED"
    df.loc[df["contam_salmo_perc"] > 0.40, "contam_species"] = "SALMON"
    df.loc[df["contam_human_perc"] > org_threshold, "contam_species"] = "HUMAN"
    df.loc[df["contam_mouse_perc"] > org_threshold, "contam_species"] = "MOUSE"

    df["contam_species_purity"] = "CLEAN"
    df.loc[df["contam_eukar_perc"] < clean_min, "contam_species_purity"] = "CONTAMINATED"
    df.loc[df["contam_eukar_perc"] < negative_min, "contam_species_purity"] = "NEGATIVE"

    # contamination flag: contaminated unless clearly HUMAN and not NEGATIVE
    df["is_contaminated"] = ~(
        (df["contam_species"] == "HUMAN") &
        (df["contam_species_purity"] != "NEGATIVE")
    )

    return df


def main():
    args = parse_args()

    # support legacy flags
    metrics_csv = args.metrics_csv or args.input
    if not metrics_csv:
        raise SystemExit("--metrics-csv (or --input) is required")
    metrics_yaml = args.metrics_yaml
    if not metrics_yaml:
        raise SystemExit("--metrics-yaml is required")

    # read CSV.gz (pandas handles gz by extension)
    df = pd.read_csv(metrics_csv)

    # compute new columns
    df = compute_derived(
        df,
        org_threshold=args.org_threshold,
        clean_min=args.euk_fraction_clean_min,
        negative_min=args.euk_fraction_negative_min,
    )

    # write CSV back in place (compressed)
    tmp_csv = metrics_csv + ".tmp"
    df.to_csv(tmp_csv, index=False, compression="gzip")
    os.replace(tmp_csv, metrics_csv)

    # update YAML schema
    schema = load_yaml(metrics_yaml)
    ensure_yaml_columns(schema, NEW_YAML_COLUMNS)
    save_yaml(schema, metrics_yaml)

    print(f"Updated CSV:  {metrics_csv}")
    print(f"Updated YAML: {metrics_yaml}")


if __name__ == "__main__":
    # deps: pandas, pyyaml
    # pip install pandas pyyaml
    main()
