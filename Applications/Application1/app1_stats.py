#!/usr/bin/env python3
# Role: Fisher's exact test comparing SAS vs EUR carrier frequency at the chr8 locus (Application 1).
import argparse
import csv
import json
import math
import os
import sys

TOOLS = ["CN.MOPS", "CNVKIT", "FREEC", "READDEPTH"]
FRAC_TO_NTOOLS = {0.25: 1, 0.5: 2, 0.75: 3, 1.0: 4}

def parse_region(s):
    chrom, rest = s.split(":")
    start, end = rest.replace(",", "").split("-")
    return chrom, int(start), int(end)

def load_populations(metadata_path):
    with open(metadata_path) as f:
        meta = json.load(f)
    pop2samples = {}
    for pop, samples in meta.items():
        pop2samples[pop] = [s for s in samples if not s.startswith("merged_")]
    return pop2samples

def find_raw_bed(data_dir, sample, tool, resolution):
    sample_dir = os.path.join(data_dir, sample)
    if tool == "READDEPTH":
        path = os.path.join(sample_dir, f"READDEPTH.{sample}.sorted.bed")
    else:
        path = os.path.join(sample_dir, f"{tool}.{resolution}.{sample}.sorted.bed")
    return path if os.path.isfile(path) else None

def find_merge_bed(merge_dir, sample, resolution, cnv_type):
    path = os.path.join(merge_dir, f"{sample}.tools_{resolution}.{cnv_type}.merge_sort.bed")
    return path if os.path.isfile(path) else None

def bed_overlaps_region(bed_path, chrom, start, end, min_support=None):
    if bed_path is None or not os.path.isfile(bed_path):
        return False
    with open(bed_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            if parts[0] != chrom:
                continue
            s, e = int(parts[1]), int(parts[2])
            if max(s, start) >= min(e, end):
                continue
            if min_support is not None:
                if len(parts) < 5:
                    continue
                try:
                    frac = round(float(parts[4]), 2)
                except ValueError:
                    continue
                level = FRAC_TO_NTOOLS.get(frac)
                if level is None or level < min_support:
                    continue
            return True
    return False

def is_carrier_tool_all(sample, tool, data_dir, resolution, chrom, start, end):
    raw = find_raw_bed(data_dir, sample, tool, resolution)
    return bed_overlaps_region(raw, chrom, start, end)

def is_carrier_cfk(sample, min_support, merge_dir, resolution, chrom, start, end):
    del_bed = find_merge_bed(merge_dir, sample, resolution, "DEL")
    dup_bed = find_merge_bed(merge_dir, sample, resolution, "DUP")
    return (bed_overlaps_region(del_bed, chrom, start, end, min_support=min_support) or
            bed_overlaps_region(dup_bed, chrom, start, end, min_support=min_support))

def fisher_exact_two_sided(a, b, c, d):
    row1, row2 = a + b, c + d
    col1, col2 = a + c, b + d
    n = row1 + row2

    def hyper_pmf(x):
        if x < 0 or x > col1 or (row1 - x) < 0 or (row1 - x) > col2:
            return 0.0
        return (math.comb(col1, x) * math.comb(col2, row1 - x)) / math.comb(n, row1)

    obs_p = hyper_pmf(a)
    x_min = max(0, row1 - col2)
    x_max = min(row1, col1)
    p_value = 0.0
    for x in range(x_min, x_max + 1):
        p_x = hyper_pmf(x)
        if p_x <= obs_p * (1 + 1e-7):
            p_value += p_x
    return min(p_value, 1.0)

def odds_ratio_ci(a, b, c, d, alpha=0.05):
    if 0 in (a, b, c, d):
        a2, b2, c2, d2 = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    else:
        a2, b2, c2, d2 = a, b, c, d
    or_val = (a2 * d2) / (b2 * c2)
    se_log_or = math.sqrt(1 / a2 + 1 / b2 + 1 / c2 + 1 / d2)
    z = 1.959963985  # 95%
    log_or = math.log(or_val)
    lo = math.exp(log_or - z * se_log_or)
    hi = math.exp(log_or + z * se_log_or)
    return or_val, lo, hi

def compute_2x2(sas_samples, eur_samples, carrier_fn):
    a = sum(1 for s in sas_samples if carrier_fn(s))          # SAS carrier
    b = len(sas_samples) - a                                   # SAS non-carrier
    c = sum(1 for s in eur_samples if carrier_fn(s))          # EUR carrier
    d = len(eur_samples) - c                                   # EUR non-carrier
    return a, b, c, d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region", required=True, help="chr8:136600000-137000000")
    ap.add_argument("--sas-populations", required=True)
    ap.add_argument("--eur-populations", required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--merge-dir", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--resolution", default="25000")
    ap.add_argument("--out-dir", default="./app1_stats_results")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    chrom, start, end = parse_region(args.region)

    pop2samples = load_populations(args.metadata)
    sas_samples = [s for p in args.sas_populations.split(",") for s in pop2samples.get(p, [])]
    eur_samples = [s for p in args.eur_populations.split(",") for s in pop2samples.get(p, [])]

    print(f"Region: {chrom}:{start}-{end}", file=sys.stderr)
    print(f"SAS pooled N={len(sas_samples)} ({args.sas_populations})", file=sys.stderr)
    print(f"EUR pooled N={len(eur_samples)} ({args.eur_populations})", file=sys.stderr)

    metrics = []
    for tool in TOOLS:
        metrics.append((
            f"{tool}-ALL",
            lambda s, t=tool: is_carrier_tool_all(s, t, args.data_dir, args.resolution, chrom, start, end)
        ))
    for k in (1, 2, 3, 4):
        metrics.append((
            f"CF>={k}",
            lambda s, k=k: is_carrier_cfk(s, k, args.merge_dir, args.resolution, chrom, start, end)
        ))

    out_path = os.path.join(args.out_dir, "app1_fisher_results.csv")
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "SAS_carrier", "SAS_total", "SAS_pct",
                          "EUR_carrier", "EUR_total", "EUR_pct",
                          "odds_ratio", "OR_95CI_low", "OR_95CI_high", "fisher_p"])

        for name, carrier_fn in metrics:
            a, b, c, d = compute_2x2(sas_samples, eur_samples, carrier_fn)
            sas_total, eur_total = a + b, c + d
            sas_pct = 100 * a / sas_total if sas_total else float("nan")
            eur_pct = 100 * c / eur_total if eur_total else float("nan")
            p = fisher_exact_two_sided(a, b, c, d)
            or_val, lo, hi = odds_ratio_ci(a, b, c, d)

            writer.writerow([name, a, sas_total, f"{sas_pct:.2f}",
                              c, eur_total, f"{eur_pct:.2f}",
                              f"{or_val:.3f}", f"{lo:.3f}", f"{hi:.3f}", f"{p:.6f}"])
            print(f"{name:12s}  SAS {a}/{sas_total} ({sas_pct:5.2f}%)  "
                  f"EUR {c}/{eur_total} ({eur_pct:5.2f}%)  "
                  f"OR={or_val:.2f} [{lo:.2f}-{hi:.2f}]  p={p:.4g}", file=sys.stderr)

    print(f"\nSaved to: {out_path}")

if __name__ == "__main__":
    main()
