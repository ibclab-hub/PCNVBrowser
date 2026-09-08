#!/usr/bin/env python3
# Role: Computes breakpoint offset between matched truth and call intervals by size bin (Supplementary Table S7).
import argparse
import csv
import os
import statistics
import subprocess
import sys
import tempfile

SIZE_BINS = [
    (50, 500, "50-500bp"),
    (500, 1000, "500-1000bp"),
    (1000, 5000, "1000-5000bp"),
    (5000, 25000, "5000-25000bp"),
    (25000, None, ">=25000bp"),
]
TOLERANCES = [100, 500, 1000, 5000]

def bin_label(size):
    for lo, hi, label in SIZE_BINS:
        if size >= lo and (hi is None or size < hi):
            return label
    return None

def find_raw_bed(data_dir, sample, tool, resolution):
    sample_dir = os.path.join(data_dir, sample)
    if tool == "READDEPTH":
        path = os.path.join(sample_dir, f"READDEPTH.{sample}.sorted.bed")
    else:
        path = os.path.join(sample_dir, f"{tool}.{resolution}.{sample}.sorted.bed")
    return path if os.path.isfile(path) else None

def extract_del_sorted(raw_path, out_path):
    rows = []
    with open(raw_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            if parts[3] != "DEL":
                continue
            rows.append((parts[0], int(parts[1]), int(parts[2])))
    rows.sort()
    with open(out_path, "w") as out:
        for chrom, start, end in rows:
            out.write(f"{chrom}\t{start}\t{end}\n")
    return len(rows)

def build_union_del(data_dir, sample, resolution, tmpdir):
    tools = ["CN.MOPS", "CNVKIT", "FREEC", "READDEPTH"]
    rows = set()
    for t in tools:
        raw = find_raw_bed(data_dir, sample, t, resolution)
        if not raw:
            continue
        with open(raw) as f:
            for line in f:
                parts = line.split()
                if len(parts) < 4 or parts[3] != "DEL":
                    continue
                rows.add((parts[0], int(parts[1]), int(parts[2])))
    out_path = os.path.join(tmpdir, f"{sample}_ALL_DEL.bed")
    with open(out_path, "w") as out:
        for chrom, start, end in sorted(rows):
            out.write(f"{chrom}\t{start}\t{end}\n")
    return out_path, len(rows)

def best_match_pairs(truth_path, call_path, reciprocal_frac):
    try:
        result = subprocess.run(
            ["bedtools", "intersect", "-a", truth_path, "-b", call_path,
             "-f", str(reciprocal_frac), "-r", "-wa", "-wb"],
            capture_output=True, text=True, check=True,
        )
    except subprocess.CalledProcessError:
        return []

    best = {}
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        t_chrom, t_start, t_end = parts[0], int(parts[1]), int(parts[2])
        c_chrom, c_start, c_end = parts[3], int(parts[4]), int(parts[5])
        start_offset = abs(t_start - c_start)
        end_offset = abs(t_end - c_end)
        total_offset = start_offset + end_offset
        key = (t_chrom, t_start, t_end)
        if key not in best or total_offset < best[key][-1]:
            best[key] = (t_chrom, t_start, t_end, c_start, c_end,
                         start_offset, end_offset, total_offset)
    return list(best.values())

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth-dir", required=True,
                     help="build_truth_beds.py / build_hprc_truth_beds.py "
                          "( DEL.bed )")
    ap.add_argument("--data-dir", required=True,
                     help="PCNVBrowser per-sample raw call bed ")
    ap.add_argument("--samples", required=True,
                     help="ID ")
    ap.add_argument("--resolution", required=True, type=int)
    ap.add_argument("--tool", required=True,
                     help="CN.MOPS / CNVKIT / FREEC / READDEPTH / ALL_4TOOL_UNION")
    ap.add_argument("--reciprocal-frac", type=float, default=0.5,
                     help="reciprocal overlap ( 0.5, "
                          "recall_by_truth_size.py )")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    samples = [s.strip() for s in args.samples.split(",") if s.strip()]

    pairs_path = os.path.join(args.out_dir, f"breakpoint_offset_pairs_{args.tool}.csv")
    summary_path = os.path.join(args.out_dir, f"breakpoint_offset_summary_{args.tool}.csv")

    all_rows = []
    with tempfile.TemporaryDirectory() as tmpdir, \
         open(pairs_path, "w", newline="") as f_pairs:

        writer = csv.writer(f_pairs)
        writer.writerow(["sample", "size_bin", "truth_start", "truth_end",
                          "call_start", "call_end", "start_offset", "end_offset",
                          "total_offset"])

        for sample in samples:
            truth_del = os.path.join(args.truth_dir, sample, "DEL.bed")
            if not os.path.isfile(truth_del):
                print(f"[SKIP] truth : {sample}", file=sys.stderr)
                continue

            if args.tool == "ALL_4TOOL_UNION":
                call_bed, n_calls = build_union_del(args.data_dir, sample,
                                                     args.resolution, tmpdir)
            else:
                raw = find_raw_bed(args.data_dir, sample, args.tool, args.resolution)
                if raw is None:
                    print(f"[SKIP] call bed : {sample} / {args.tool}", file=sys.stderr)
                    continue
                call_bed = os.path.join(tmpdir, f"{sample}_{args.tool}_DEL.bed")
                extract_del_sorted(raw, call_bed)

            matched = best_match_pairs(truth_del, call_bed, args.reciprocal_frac)
            for (chrom, t_start, t_end, c_start, c_end,
                 start_off, end_off, total_off) in matched:
                size = t_end - t_start
                label = bin_label(size)
                if label is None:
                    continue
                writer.writerow([sample, label, t_start, t_end, c_start, c_end,
                                  start_off, end_off, total_off])
                all_rows.append((label, total_off))

            print(f"{sample} : matched pairs={len(matched)}", file=sys.stderr)

    # size_bin 
    from collections import defaultdict
    by_bin = defaultdict(list)
    for label, total_off in all_rows:
        by_bin[label].append(total_off)

    with open(summary_path, "w", newline="") as f_sum:
        writer = csv.writer(f_sum)
        header = ["size_bin", "n_pairs", "median_total_offset_bp", "mean_total_offset_bp"]
        header += [f"pct_within_{t}bp" for t in TOLERANCES]
        writer.writerow(header)
        for _, _, label in SIZE_BINS:
            offs = by_bin.get(label, [])
            if not offs:
                writer.writerow([label, 0, "NA", "NA"] + ["NA"] * len(TOLERANCES))
                continue
            row = [label, len(offs),
                   f"{statistics.median(offs):.1f}",
                   f"{statistics.mean(offs):.1f}"]
            for t in TOLERANCES:
                pct = 100.0 * sum(1 for o in offs if o <= t) / len(offs)
                row.append(f"{pct:.1f}")
            writer.writerow(row)

    print(f"(A) pair-level: {pairs_path}", file=sys.stderr)
    print(f"(B) size-bin summary: {summary_path}", file=sys.stderr)

if __name__ == "__main__":
    main()
