#!/usr/bin/env python3
# Role: Computes pairwise Jaccard index and reciprocal overlap between the four CNV callers (Supplementary Table S3).

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
from itertools import combinations

TOOLS = ["CN.MOPS", "CNVKIT", "FREEC", "READDEPTH"]
TYPES = ["DEL", "DUP"]

def find_bed(sample_dir, sample, tool, resolution):
    if tool == "READDEPTH":
        path = os.path.join(sample_dir, f"READDEPTH.{sample}.sorted.bed")
    else:
        path = os.path.join(sample_dir, f"{tool}.{resolution}.{sample}.sorted.bed")
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

def bedtools_jaccard(path_a, path_b):
    try:
        result = subprocess.run(
            ["bedtools", "jaccard", "-a", path_a, "-b", path_b],
            capture_output=True, text=True, check=True,
        )
    except subprocess.CalledProcessError:
        return None
    lines = result.stdout.strip().splitlines()
    if len(lines) < 2:
        return None
    fields = lines[1].split("\t")
    try:
        return float(fields[2])
    except (IndexError, ValueError):
        return None

def bedtools_reciprocal_frac(path_a, path_b, n_a, frac):
    if n_a == 0:
        return None
    try:
        result = subprocess.run(
            ["bedtools", "intersect", "-u", "-f", str(frac), "-r",
             "-a", path_a, "-b", path_b],
            capture_output=True, text=True, check=True,
        )
    except subprocess.CalledProcessError:
        return None
    n_hit = len([l for l in result.stdout.splitlines() if l.strip()])
    return n_hit / n_a

def load_population_map(metadata_path):
    with open(metadata_path) as f:
        meta = json.load(f)
    sample2pop = {}
    for pop, samples in meta.items():
        for s in samples:
            if s.startswith("merged_"):
                continue
            sample2pop[s] = pop
    return sample2pop

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True, help="")
    ap.add_argument("--sample-list", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--resolution", default="25000",
                     help="CN.MOPS/CNVKIT/FREEC bin size (5000/25000/50000). "
                          "READDEPTH .")
    ap.add_argument("--reciprocal-frac", type=float, default=0.5)
    ap.add_argument("--out-dir", default="./concordance_results")
    ap.add_argument("--limit", type=int, default=None,
                     help="N ")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    with open(args.sample_list) as f:
        samples = [l.strip() for l in f if l.strip()]
    if args.limit:
        samples = samples[: args.limit]

    sample2pop = load_population_map(args.metadata)
    pairs = list(combinations(TOOLS, 2))

    out_path = os.path.join(
        args.out_dir, f"pairwise_concordance_res{args.resolution}.csv"
    )
    missing_log_path = os.path.join(
        args.out_dir, f"missing_files_res{args.resolution}.log"
    )

    n_total = len(samples)
    with open(out_path, "w", newline="") as out_f, \
         open(missing_log_path, "w") as miss_f, \
         tempfile.TemporaryDirectory() as tmpdir:

        writer = csv.writer(out_f)
        writer.writerow([
            "sample", "population", "cnv_type", "tool_a", "tool_b",
            "n_calls_a", "n_calls_b", "jaccard",
            f"recip{int(args.reciprocal_frac*100)}_a_in_b",
            f"recip{int(args.reciprocal_frac*100)}_b_in_a",
        ])

        for idx, sample in enumerate(samples, 1):
            sample_dir = os.path.join(args.data_dir, sample)
            pop = sample2pop.get(sample, "NA")

            for cnv_type in TYPES:
                tmp_paths = {}
                n_calls = {}
                for tool in TOOLS:
                    src = find_bed(sample_dir, sample, tool, args.resolution)
                    if src is None:
                        miss_f.write(f"{sample}\t{tool}\t{cnv_type}\tFILE_NOT_FOUND\n")
                        tmp_paths[tool] = None
                        n_calls[tool] = 0
                        continue
                    dst = os.path.join(tmpdir, f"{tool}_{cnv_type}.bed")
                    n = extract_type_sorted(src, cnv_type, dst)
                    tmp_paths[tool] = dst if n > 0 else None
                    n_calls[tool] = n

                for tool_a, tool_b in pairs:
                    pa, pb = tmp_paths[tool_a], tmp_paths[tool_b]
                    na, nb = n_calls[tool_a], n_calls[tool_b]
                    if pa is None and pb is None:
                        continue # , 
                    if pa is None or pb is None:
                        j = 0.0
                        ra = 0.0 if na > 0 else None
                        rb = 0.0 if nb > 0 else None
                    else:
                        j = bedtools_jaccard(pa, pb)
                        ra = bedtools_reciprocal_frac(pa, pb, na, args.reciprocal_frac)
                        rb = bedtools_reciprocal_frac(pb, pa, nb, args.reciprocal_frac)

                    writer.writerow([
                        sample, pop, cnv_type, tool_a, tool_b, na, nb,
                        f"{j:.4f}" if j is not None else "NA",
                        f"{ra:.4f}" if ra is not None else "NA",
                        f"{rb:.4f}" if rb is not None else "NA",
                    ])

            if idx % 100 == 0 or idx == n_total:
                print(f"[{idx}/{n_total}] {sample} ", file=sys.stderr)

    print(f". : {out_path}")
    print(f" : {missing_log_path}")

if __name__ == "__main__":
    main()
