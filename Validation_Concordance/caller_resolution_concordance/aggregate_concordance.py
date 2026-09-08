#!/usr/bin/env python3
# Role: Aggregates per-sample pairwise caller concordance into summary statistics (Supplementary Table S3).
import argparse
import os
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    for col in ["jaccard"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    recip_cols = [c for c in df.columns if c.startswith("recip")]
    for col in recip_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    # recip_a_in_b / recip_b_in_a symmetric 
    df["recip_mean"] = df[recip_cols].mean(axis=1, skipna=True)

    os.makedirs(args.out_dir, exist_ok=True)

    # 1) : tool-pair x cnv_type
    overall = (
        df.groupby(["cnv_type", "tool_a", "tool_b"])
        .agg(
            n_samples=("sample", "nunique"),
            jaccard_mean=("jaccard", "mean"),
            jaccard_median=("jaccard", "median"),
            jaccard_std=("jaccard", "std"),
            recip_mean=("recip_mean", "mean"),
            recip_median=("recip_mean", "median"),
        )
        .reset_index()
        .sort_values(["cnv_type", "jaccard_mean"], ascending=[True, False])
    )
    overall_path = os.path.join(args.out_dir, "summary_overall.csv")
    overall.to_csv(overall_path, index=False)

    # 2) population: tool-pair x cnv_type x population
    by_pop = (
        df.groupby(["cnv_type", "tool_a", "tool_b", "population"])
        .agg(
            n_samples=("sample", "nunique"),
            jaccard_mean=("jaccard", "mean"),
            recip_mean=("recip_mean", "mean"),
        )
        .reset_index()
        .sort_values(["cnv_type", "tool_a", "tool_b", "population"])
    )
    by_pop_path = os.path.join(args.out_dir, "summary_by_population.csv")
    by_pop.to_csv(by_pop_path, index=False)

    print("=== () ===")
    print(overall.to_string(index=False))
    print(f"\nSaved to: {overall_path}")
    print(f"population : {by_pop_path}")

if __name__ == "__main__":
    main()
