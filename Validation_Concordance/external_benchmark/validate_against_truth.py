#!/usr/bin/env python3
# Role: Computes recall and precision of caller-specific and support-level-filtered CNV calls against truth sets (Supplementary Tables S5, S6).
import argparse
import csv
import os
import subprocess
import sys
import tempfile

TOOLS = ["CN.MOPS", "CNVKIT", "FREEC", "READDEPTH"]
FRAC_TO_NTOOLS = {0.25: 1, 0.5: 2, 0.75: 3, 1.0: 4}

def find_raw_bed(sample_dir, sample, tool, resolution):
    if tool == "READDEPTH":
        path = os.path.join(sample_dir, f"READDEPTH.{sample}.sorted.bed")
    else:
        path = os.path.join(sample_dir, f"{tool}.{resolution}.{sample}.sorted.bed")
    return path if os.path.isfile(path) else None

def find_merge_bed(merge_dir, sample, resolution, cnv_type):
    path = os.path.join(merge_dir, f"{sample}.tools_{resolution}.{cnv_type}.merge_sort.bed")
    return path if os.path.isfile(path) else None

def extract_type_sorted(src_path, cnv_type, dst_path):
    n = 0
    rows = []
    with open(src_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            chrom, start, end, typ = parts[0], parts[1], parts[2], parts[3]
            if typ != cnv_type:
                continue
            rows.append((chrom, int(start), int(end)))
            n += 1
    rows.sort()
    with open(dst_path, "w") as out:
        for chrom, start, end in rows:
            out.write(f"{chrom}\t{start}\t{end}\n")
    return n

def merge_by_support_threshold(merge_path, min_ntools, dst_path):
    rows = []
    with open(merge_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 5:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            try:
                frac = round(float(parts[4]), 2)
            except ValueError:
                continue
            level = FRAC_TO_NTOOLS.get(frac)
            if level is not None and level >= min_ntools:
                rows.append((chrom, start, end))
    rows.sort()
    n = 0
    with open(dst_path, "w") as out:
        for chrom, start, end in rows:
            out.write(f"{chrom}\t{start}\t{end}\n")
            n += 1
    return n

def restrict_to_confident(call_bed, confident_bed, dst_path):
    if confident_bed is None or not os.path.isfile(confident_bed):
        return call_bed
    try:
        result = subprocess.run(
            ["bedtools", "intersect", "-u", "-a", call_bed, "-b", confident_bed],
            capture_output=True, text=True, check=True,
        )
    except subprocess.CalledProcessError:
        return call_bed
    lines = [l for l in result.stdout.splitlines() if l.strip()]
    if not lines:
        return None
    with open(dst_path, "w") as f:
        for l in lines:
            f.write(l + "\n")
    return dst_path

def find_confident_bed(truth_dir, sample):
    path = os.path.join(truth_dir, sample, "CONFIDENT.bed")
    return path if os.path.isfile(path) else None

def sort_bed(path):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return path
    tmp = path + ".sorted"
    subprocess.run(f"sort -k1,1 -k2,2n {path} > {tmp}", shell=True, check=True)
    os.replace(tmp, path)
    return path

def n_lines(path):
    if not os.path.isfile(path):
        return 0
    with open(path) as f:
        return sum(1 for _ in f)

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

def combine_beds(paths, dst_path):
    rows = []
    for p in paths:
        if p and os.path.isfile(p):
            with open(p) as f:
                for line in f:
                    if line.strip():
                        rows.append(line.strip())
    rows = sorted(set(rows))
    with open(dst_path, "w") as out:
        for r in rows:
            out.write(r + "\n")
    return len(rows)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth-dir", required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--merge-dir", required=True)
    ap.add_argument("--samples", required=True, help="ID ")
    ap.add_argument("--resolution", default="25000")
    ap.add_argument("--reciprocal-frac", type=float, default=0.5)
    ap.add_argument("--include-ins", action="store_true",
                     help="DUP truth INS.bed ")
    ap.add_argument("--out-dir", default="./validation_results")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    samples = [s.strip() for s in args.samples.split(",") if s.strip()]

    per_tool_path = os.path.join(args.out_dir, f"recall_precision_per_tool_res{args.resolution}.csv")
    per_support_path = os.path.join(args.out_dir, f"recall_precision_by_support_res{args.resolution}.csv")
    skipped_log = os.path.join(args.out_dir, "skipped_samples.log")

    with open(per_tool_path, "w", newline="") as f1, \
         open(per_support_path, "w", newline="") as f2, \
         open(skipped_log, "w") as skip_f, \
         tempfile.TemporaryDirectory() as tmpdir:

        w1 = csv.writer(f1)
        w1.writerow(["sample", "tool", "cnv_type", "n_calls", "n_truth", "recall", "precision"])
        w2 = csv.writer(f2)
        w2.writerow(["sample", "cnv_type", "min_ntools_support", "n_calls", "n_truth", "recall", "precision"])

        for sample in samples:
            truth_del = os.path.join(args.truth_dir, sample, "DEL.bed")
            truth_dup_parts = []
            truth_dup_direct = os.path.join(args.truth_dir, sample, "DUP.bed")
            if os.path.isfile(truth_dup_direct):
                truth_dup_parts.append(truth_dup_direct)
            if args.include_ins:
                truth_ins = os.path.join(args.truth_dir, sample, "INS.bed")
                if os.path.isfile(truth_ins):
                    truth_dup_parts.append(truth_ins)

            truth_dup = None
            if truth_dup_parts:
                truth_dup = os.path.join(tmpdir, f"{sample}_truth_DUP.bed")
                combine_beds(truth_dup_parts, truth_dup)
                sort_bed(truth_dup)

            if os.path.isfile(truth_del):
                sort_bed(truth_del)
            else:
                truth_del = None

            if truth_del is None and truth_dup is None:
                skip_f.write(f"{sample}\tNO_TRUTH_FOUND\n")
                continue

            sample_dir = os.path.join(args.data_dir, sample)
            truth_by_type = {"DEL": truth_del, "DUP": truth_dup}

            # (A) per-tool
            for tool in TOOLS:
                src = find_raw_bed(sample_dir, sample, tool, args.resolution)
                if src is None:
                    skip_f.write(f"{sample}\t{tool}\tRAW_BED_NOT_FOUND\n")
                    continue
                for cnv_type in ("DEL", "DUP"):
                    truth_bed = truth_by_type[cnv_type]
                    if truth_bed is None:
                        continue
                    call_bed = os.path.join(tmpdir, f"{sample}_{tool}_{cnv_type}.bed")
                    n_calls = extract_type_sorted(src, cnv_type, call_bed)
                    confident_bed = find_confident_bed(args.truth_dir, sample)
                    if confident_bed and n_calls > 0:
                        restricted = restrict_to_confident(
                            call_bed, confident_bed,
                            os.path.join(tmpdir, f"{sample}_{tool}_{cnv_type}_conf.bed")
                        )
                        if restricted is None:
                            n_calls = 0
                        else:
                            call_bed = restricted
                            n_calls = n_lines(call_bed)
                    n_truth = n_lines(truth_bed)
                    if n_calls == 0 or n_truth == 0:
                        recall = precision = None
                    else:
                        recall = reciprocal_frac(truth_bed, call_bed, n_truth, args.reciprocal_frac)
                        precision = reciprocal_frac(call_bed, truth_bed, n_calls, args.reciprocal_frac)
                    w1.writerow([sample, tool, cnv_type, n_calls, n_truth,
                                 f"{recall:.4f}" if recall is not None else "NA",
                                 f"{precision:.4f}" if precision is not None else "NA"])

            # (B) support-level threshold
            for cnv_type in ("DEL", "DUP"):
                truth_bed = truth_by_type[cnv_type]
                if truth_bed is None:
                    continue
                merge_path = find_merge_bed(args.merge_dir, sample, args.resolution, cnv_type)
                if merge_path is None:
                    skip_f.write(f"{sample}\tMERGE\t{cnv_type}\tMERGE_FILE_NOT_FOUND\n")
                    continue
                n_truth = n_lines(truth_bed)
                for min_n in (1, 2, 3, 4):
                    call_bed = os.path.join(tmpdir, f"{sample}_{cnv_type}_min{min_n}.bed")
                    n_calls = merge_by_support_threshold(merge_path, min_n, call_bed)
                    confident_bed = find_confident_bed(args.truth_dir, sample)
                    if confident_bed and n_calls > 0:
                        restricted = restrict_to_confident(
                            call_bed, confident_bed,
                            os.path.join(tmpdir, f"{sample}_{cnv_type}_min{min_n}_conf.bed")
                        )
                        if restricted is None:
                            n_calls = 0
                        else:
                            call_bed = restricted
                            n_calls = n_lines(call_bed)
                    if n_calls == 0 or n_truth == 0:
                        recall = precision = None
                    else:
                        recall = reciprocal_frac(truth_bed, call_bed, n_truth, args.reciprocal_frac)
                        precision = reciprocal_frac(call_bed, truth_bed, n_calls, args.reciprocal_frac)
                    w2.writerow([sample, cnv_type, min_n, n_calls, n_truth,
                                 f"{recall:.4f}" if recall is not None else "NA",
                                 f"{precision:.4f}" if precision is not None else "NA"])

            print(f"{sample} ", file=sys.stderr)

    print(f"(A) tool : {per_tool_path}")
    print(f"(B) support-level : {per_support_path}")
    print(f" : {skipped_log}")

if __name__ == "__main__":
    main()
