#!/usr/bin/env python3
# Role: Computes size-stratified recall/precision by truth-call size bin (Supplementary Table S6).
import argparse
import csv
import os
import subprocess
import sys
import tempfile

SIZE_BINS = [(50, 500), (500, 1000), (1000, 5000), (5000, 25000), (25000, float("inf"))]

def bin_label(lo, hi):
    if hi == float("inf"):
        return f">={lo}bp"
    return f"{lo}-{hi}bp"

def find_raw_bed(data_dir, sample, tool, resolution):
    sample_dir = os.path.join(data_dir, sample)
    if tool == "READDEPTH":
        path = os.path.join(sample_dir, f"READDEPTH.{sample}.sorted.bed")
    else:
        path = os.path.join(sample_dir, f"{tool}.{resolution}.{sample}.sorted.bed")
    return path if os.path.isfile(path) else None

def extract_del_sorted(src_path, dst_path):
    rows = []
    with open(src_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 4 or parts[3] != "DEL":
                continue
            rows.append((parts[0], int(parts[1]), int(parts[2])))
    rows.sort()
    with open(dst_path, "w") as out:
        for chrom, start, end in rows:
            out.write(f"{chrom}\t{start}\t{end}\n")
    return len(rows)

def split_truth_by_size(truth_del_path, dst_dir, sample):
    bins = {b: [] for b in SIZE_BINS}
    with open(truth_del_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            start, end = int(parts[1]), int(parts[2])
            size = end - start
            for lo, hi in SIZE_BINS:
                if lo <= size < hi:
                    bins[(lo, hi)].append((parts[0], start, end))
                    break
    out = {}
    for (lo, hi), rows in bins.items():
        rows.sort()
        path = os.path.join(dst_dir, f"{sample}_truth_{lo}_{hi}.bed")
        with open(path, "w") as f:
            for chrom, s, e in rows:
                f.write(f"{chrom}\t{s}\t{e}\n")
        out[(lo, hi)] = (path, len(rows))
    return out

def split_calls_by_size(call_bed_path, dst_dir, sample, tag):
    bins = {b: [] for b in SIZE_BINS}
    with open(call_bed_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            start, end = int(parts[1]), int(parts[2])
            size = end - start
            for lo, hi in SIZE_BINS:
                if lo <= size < hi:
                    bins[(lo, hi)].append((parts[0], start, end))
                    break
    out = {}
    for (lo, hi), rows in bins.items():
        rows.sort()
        path = os.path.join(dst_dir, f"{sample}_{tag}_calls_{lo}_{hi}.bed")
        with open(path, "w") as f:
            for chrom, s, e in rows:
                f.write(f"{chrom}\t{s}\t{e}\n")
        out[(lo, hi)] = (path, len(rows))
    return out

def reciprocal_frac(path_a, path_b, n_a, frac=0.5):
    if n_a == 0:
        return None
    try:
        result = subprocess.run(
            ["bedtools", "intersect", "-u", "-f", str(frac), "-r", "-a", path_a, "-b", path_b],
            capture_output=True, text=True, check=True,
        )
    except subprocess.CalledProcessError:
        return None
    n_hit = len([l for l in result.stdout.splitlines() if l.strip()])
    return n_hit / n_a

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth-dir", required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--resolution", default="25000")
    ap.add_argument("--tool", default="READDEPTH",
                     help="tool raw call recall ( READDEPTH, "
                          " tool 'ALL_4TOOL_UNION' )")
    ap.add_argument("--reciprocal-frac", type=float, default=0.5)
    ap.add_argument("--out-dir", default="./size_stratified_results")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    samples = [s.strip() for s in args.samples.split(",") if s.strip()]

    out_path = os.path.join(args.out_dir, f"recall_by_size_{args.tool}.csv")
    precision_out_path = os.path.join(args.out_dir, f"precision_by_size_{args.tool}.csv")

    with open(out_path, "w", newline="") as f, \
         open(precision_out_path, "w", newline="") as fp, \
         tempfile.TemporaryDirectory() as tmpdir:

        writer = csv.writer(f)
        writer.writerow(["sample", "size_bin", "n_truth_in_bin", "n_calls", "recall"])
        writer_p = csv.writer(fp)
        writer_p.writerow(["sample", "size_bin", "n_calls_in_bin", "n_truth_total", "precision"])

        for sample in samples:
            truth_del = os.path.join(args.truth_dir, sample, "DEL.bed")
            if not os.path.isfile(truth_del):
                continue

            size_bins = split_truth_by_size(truth_del, tmpdir, sample)
            n_truth_total = sum(n for _, n in size_bins.values())

            if args.tool == "ALL_4TOOL_UNION":
                tools = ["CN.MOPS", "CNVKIT", "FREEC", "READDEPTH"]
                paths = []
                for t in tools:
                    raw = find_raw_bed(args.data_dir, sample, t, args.resolution)
                    if raw:
                        p = os.path.join(tmpdir, f"{sample}_{t}_DEL.bed")
                        extract_del_sorted(raw, p)
                        paths.append(p)
                combined = os.path.join(tmpdir, f"{sample}_ALL_DEL.bed")
                rows = set()
                for p in paths:
                    with open(p) as pf:
                        for l in pf:
                            if l.strip():
                                rows.add(l.strip())
                with open(combined, "w") as cf:
                    for r in sorted(rows):
                        cf.write(r + "\n")
                call_bed = combined
                n_calls = len(rows)
            else:
                raw = find_raw_bed(args.data_dir, sample, args.tool, args.resolution)
                if raw is None:
                    continue
                call_bed = os.path.join(tmpdir, f"{sample}_{args.tool}_DEL.bed")
                n_calls = extract_del_sorted(raw, call_bed)

            # (A) recall: truth , (call_bed) 
            for (lo, hi), (truth_bin_path, n_truth_bin) in size_bins.items():
                label = bin_label(lo, hi)
                if n_truth_bin == 0 or n_calls == 0:
                    writer.writerow([sample, label, n_truth_bin, n_calls, "NA"])
                    continue
                recall = reciprocal_frac(truth_bin_path, call_bed, n_truth_bin, args.reciprocal_frac)
                writer.writerow([sample, label, n_truth_bin,
                                  n_calls, f"{recall:.4f}" if recall is not None else "NA"])

            # (B) precision: , truth (truth_del) 
            if n_calls > 0:
                call_bins = split_calls_by_size(call_bed, tmpdir, sample, args.tool)
                for (lo, hi), (call_bin_path, n_calls_bin) in call_bins.items():
                    label = bin_label(lo, hi)
                    if n_calls_bin == 0 or n_truth_total == 0:
                        writer_p.writerow([sample, label, n_calls_bin, n_truth_total, "NA"])
                        continue
                    precision = reciprocal_frac(call_bin_path, truth_del, n_calls_bin,
                                                 args.reciprocal_frac)
                    writer_p.writerow([sample, label, n_calls_bin, n_truth_total,
                                        f"{precision:.4f}" if precision is not None else "NA"])

            print(f"{sample} ", file=sys.stderr)

    print(f"\nrecall : {out_path}")
    print(f"precision : {precision_out_path}")

if __name__ == "__main__":
    main()
