#!/usr/bin/env python3
"""
aggregate_recall_precision.py

recall_by_truth_size_fixed.py가 tool별로 만든
recall_by_size_{TOOL}.csv / precision_by_size_{TOOL}.csv (샘플별 raw 값)를
size bin별로 macro-average(샘플 평균)해서 논문에 들어갈 요약 표를 만든다.

- recall: sample당 recall 값이 있는 행만 평균 (NA 제외)
- precision: sample당 precision 값이 있는 행만 평균 (NA 제외)
- ">=25000bp" bin이 본문 텍스트("recall for deletions >=25 kb ...")에 해당하는 행

사용 예:
  python3 aggregate_recall_precision.py \
      --results-dir ./size_stratified_results_25k_confident \
      --tools READDEPTH,CN.MOPS,CNVKIT,FREEC,ALL_4TOOL_UNION \
      --out summary_25k_confident.csv
"""
import argparse
import csv
import os
import statistics


def load_metric_csv(path, metric_col):
    """size_bin -> list of (sample, value) 반환, NA/빈 값은 제외."""
    by_bin = {}
    if not os.path.isfile(path):
        return by_bin
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            size_bin = row["size_bin"]
            val = row[metric_col]
            if val is None or val == "" or val == "NA":
                continue
            try:
                v = float(val)
            except ValueError:
                continue
            by_bin.setdefault(size_bin, []).append((row["sample"], v))
    return by_bin


def summarize(values):
    if not values:
        return None, None, 0
    nums = [v for _, v in values]
    mean = statistics.mean(nums)
    return mean, (min(nums), max(nums)), len(nums)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--tools", required=True, help="콤마로 구분된 tool 이름 목록")
    ap.add_argument("--out", default="summary.csv")
    args = ap.parse_args()

    tools = [t.strip() for t in args.tools.split(",") if t.strip()]

    rows_out = []
    for tool in tools:
        recall_path = os.path.join(args.results_dir, f"recall_by_size_{tool}.csv")
        precision_path = os.path.join(args.results_dir, f"precision_by_size_{tool}.csv")

        recall_by_bin = load_metric_csv(recall_path, "recall")
        precision_by_bin = load_metric_csv(precision_path, "precision")

        all_bins = sorted(set(recall_by_bin) | set(precision_by_bin))
        for size_bin in all_bins:
            r_mean, r_range, r_n = summarize(recall_by_bin.get(size_bin, []))
            p_mean, p_range, p_n = summarize(precision_by_bin.get(size_bin, []))
            rows_out.append({
                "tool": tool,
                "size_bin": size_bin,
                "mean_recall": f"{r_mean*100:.1f}%" if r_mean is not None else "NA",
                "n_samples_recall": r_n,
                "mean_precision": f"{p_mean*100:.1f}%" if p_mean is not None else "NA",
                "n_samples_precision": p_n,
            })

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "tool", "size_bin", "mean_recall", "n_samples_recall",
            "mean_precision", "n_samples_precision"
        ])
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"저장됨: {args.out}")
    print()
    print(f"{'tool':<18}{'size_bin':<12}{'recall':<10}{'n':<5}{'precision':<10}{'n'}")
    for r in rows_out:
        if r["size_bin"] == ">=25000bp":
            print(f"{r['tool']:<18}{r['size_bin']:<12}{r['mean_recall']:<10}"
                  f"{r['n_samples_recall']:<5}{r['mean_precision']:<10}{r['n_samples_precision']}")


if __name__ == "__main__":
    main()
