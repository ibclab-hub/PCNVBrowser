#!/usr/bin/env python3
# Role: Compares alternative carrier-frequency metric definitions for internal consistency.
import argparse
import csv
import os
import sys
from collections import defaultdict

def load_bed_with_freq(path):
    out = defaultdict(list)
    if not os.path.isfile(path):
        return out
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                start, end, freq = int(parts[1]), int(parts[2]), float(parts[4])
            except ValueError:
                continue
            out[parts[0]].append((start, end, parts[3], freq))
    for c in out:
        out[c].sort()
    return out

def weighted_mean_freq(bed_dict):
    total_bp = 0
    weighted_sum = 0.0
    for chrom, rows in bed_dict.items():
        for start, end, _type, freq in rows:
            length = end - start
            total_bp += length
            weighted_sum += freq * length
    if total_bp == 0:
        return None
    return weighted_sum / total_bp

def pair_overlaps(bed_a, bed_b):
    pairs = []
    for chrom, a_rows in bed_a.items():
        b_rows = bed_b.get(chrom, [])
        for a_s, a_e, _at, a_freq in a_rows:
            for b_s, b_e, _bt, b_freq in b_rows:
                lo, hi = max(a_s, b_s), min(a_e, b_e)
                if lo < hi:
                    pairs.append((hi - lo, a_freq, b_freq))
    return pairs

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old-all-dir", required=True,
                     help="{pop}.tools_{res}.{TYPE}.merge_sort.bed (all_generator.pl ) ")
    ap.add_argument("--carrier-freq-dir", required=True,
                     help="carrier_frequency.py ")
    ap.add_argument("--resolution", default="25000")
    ap.add_argument("--populations", required=True, help="population ")
    ap.add_argument("--out-dir", default="./freq_comparison")
    ap.add_argument("--divergence-threshold", type=float, default=0.3,
                     help="old-new ' ' ")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    pops = args.populations.split(",")

    # (A) Monotonicity check
    mono_path = os.path.join(args.out_dir, "monotonicity_check.csv")
    with open(mono_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["population", "cnv_type", "min1_freq", "min2_freq",
                          "min3_freq", "min4_freq", "is_monotonic_nonincreasing"])
        for pop in pops:
            for cnv_type in ("DEL", "DUP", "ALL"):
                freqs = []
                for level in (1, 2, 3, 4):
                    path = os.path.join(
                        args.carrier_freq_dir,
                        f"{pop}.tools_{args.resolution}.{cnv_type}.carrier_freq_min{level}.bed"
                    )
                    bed = load_bed_with_freq(path)
                    freqs.append(weighted_mean_freq(bed))
                valid = [x for x in freqs if x is not None]
                is_mono = all(valid[i] >= valid[i + 1] - 1e-9 for i in range(len(valid) - 1)) if len(valid) > 1 else None
                writer.writerow([pop, cnv_type] + [f"{x:.5f}" if x is not None else "NA" for x in freqs] +
                                 [is_mono])

    # (B) ALL vs CF>=1 (type : old DEL vs new DEL CF1, old DUP vs new DUP CF1,
    #     old ALL(mixed) vs new ALL(union) CF1)
    detail_path = os.path.join(args.out_dir, "old_vs_new_divergence_detail.csv")
    summary_path = os.path.join(args.out_dir, "old_vs_new_divergence_summary.csv")

    summary_rows = []
    with open(detail_path, "w", newline="") as detail_f:
        writer = csv.writer(detail_f)
        writer.writerow(["population", "cnv_type", "overlap_bp", "old_freq", "new_cf1_freq",
                          "delta_new_minus_old"])

        for pop in pops:
            for cnv_type in ("DEL", "DUP", "ALL"):
                old_path = os.path.join(
                    args.old_all_dir, f"{pop}.tools_{args.resolution}.{cnv_type}.merge_sort.bed"
                )
                new_path = os.path.join(
                    args.carrier_freq_dir,
                    f"{pop}.tools_{args.resolution}.{cnv_type}.carrier_freq_min1.bed"
                )
                old_bed = load_bed_with_freq(old_path)
                new_bed = load_bed_with_freq(new_path)
                pairs = pair_overlaps(old_bed, new_bed)
                if not pairs:
                    continue

                deltas = []
                large_divergence_bp = 0
                total_bp = 0
                for bp, old_f, new_f in pairs:
                    delta = new_f - old_f
                    deltas.append((bp, delta))
                    total_bp += bp
                    if abs(delta) >= args.divergence_threshold:
                        large_divergence_bp += bp
                    writer.writerow([pop, cnv_type, bp, f"{old_f:.5f}", f"{new_f:.5f}", f"{delta:.5f}"])

                weighted_mean_delta = sum(bp * d for bp, d in deltas) / total_bp if total_bp else None
                pct_large_divergence = (large_divergence_bp / total_bp * 100) if total_bp else None

                summary_rows.append([
                    pop, cnv_type, total_bp,
                    f"{weighted_mean_delta:.5f}" if weighted_mean_delta is not None else "NA",
                    f"{pct_large_divergence:.2f}" if pct_large_divergence is not None else "NA",
                ])

    with open(summary_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["population", "cnv_type", "total_overlap_bp",
                          "weighted_mean_delta_new_minus_old",
                          f"pct_bp_with_abs_delta_ge_{args.divergence_threshold}"])
        writer.writerows(summary_rows)

    print(f"(A) Monotonicity check: {mono_path}")
    print(f"(B) delta: {detail_path}")
    print(f"(B) : {summary_path}")

    print("\n=== (A) Monotonicity ===")
    with open(mono_path) as f:
        print(f.read())
    print("=== (B) Old vs New divergence ===")
    with open(summary_path) as f:
        print(f.read())

if __name__ == "__main__":
    main()
