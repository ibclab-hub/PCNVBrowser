#!/usr/bin/env python3
# Role: Aggregates within-caller resolution concordance results (Supplementary Table S4).
import argparse
import os
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    df = pd.read_csv(args.input)
    df["jaccard"] = pd.to_numeric(df["jaccard"], errors="coerce")
    recip_cols = [c for c in df.columns if c.startswith("recip")]
    for c in recip_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["recip_mean"] = df[recip_cols].mean(axis=1, skipna=True)

    # (A) pairwise concordance 
    overall = (
        df.groupby(["tool", "cnv_type", "win_a", "win_b"])
        .agg(
            n_samples=("sample", "nunique"),
            jaccard_mean=("jaccard", "mean"),
            jaccard_median=("jaccard", "median"),
            jaccard_std=("jaccard", "std"),
            recip_mean=("recip_mean", "mean"),
            recip_median=("recip_mean", "median"),
        )
        .reset_index()
        .sort_values(["tool", "cnv_type", "win_a", "win_b"])
    )
    overall_path = os.path.join(args.out_dir, "summary_resolution_concordance.csv")
    overall.to_csv(overall_path, index=False)

    # (B) tool x resolution x cnv_type 
    long_a = df[["sample", "tool", "cnv_type", "win_a", "n_calls_a"]].rename(
        columns={"win_a": "win", "n_calls_a": "n_calls"}
    )
    long_b = df[["sample", "tool", "cnv_type", "win_b", "n_calls_b"]].rename(
        columns={"win_b": "win", "n_calls_b": "n_calls"}
    )
    long_all = pd.concat([long_a, long_b], ignore_index=True).drop_duplicates()

    call_counts = (
        long_all.groupby(["tool", "cnv_type", "win"])["n_calls"]
        .agg(mean="mean", median="median", std="std", n_samples="count")
        .reset_index()
        .sort_values(["tool", "cnv_type", "win"])
    )
    call_counts_path = os.path.join(args.out_dir, "summary_resolution_call_counts.csv")
    call_counts.to_csv(call_counts_path, index=False)

    print("=== (A) Resolution pairwise concordance ===")
    print(overall.to_string(index=False))
    print(f"\nSaved to: {overall_path}")

    print("\n=== (B) Tool x Resolution (DEL/DUP ) ===")
    print(call_counts.to_string(index=False))
    print(f"\nSaved to: {call_counts_path}")

if __name__ == "__main__":
    main()
