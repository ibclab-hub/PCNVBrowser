#!/usr/bin/env python3
# Role: QC check for individuals classified as carriers of both a DEL and a DUP at the same locus and support level.
import argparse
import csv
import os
import sys
from collections import Counter

TOOLS = ["CN.MOPS", "CNVKIT", "FREEC", "READDEPTH"]
FRAC_TO_NTOOLS = {0.25: 1, 0.5: 2, 0.75: 3, 1.0: 4}

def find_merge_bed(merge_dir, sample, resolution, cnv_type):
    path = os.path.join(merge_dir, f"{sample}.tools_{resolution}.{cnv_type}.merge_sort.bed")
    return path if os.path.isfile(path) else None

def find_raw_bed(data_dir, sample, tool, resolution):
    sample_dir = os.path.join(data_dir, sample)
    if tool == "READDEPTH":
        path = os.path.join(sample_dir, f"READDEPTH.{sample}.sorted.bed")
    else:
        path = os.path.join(sample_dir, f"{tool}.{resolution}.{sample}.sorted.bed")
    return path if os.path.isfile(path) else None

def load_presence_intervals(merge_bed_path, min_support):
    from collections import defaultdict
    out = defaultdict(list)
    with open(merge_bed_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                frac = round(float(parts[4]), 2)
            except ValueError:
                continue
            level = FRAC_TO_NTOOLS.get(frac)
            if level is None or level < min_support:
                continue
            out[parts[0]].append((int(parts[1]), int(parts[2])))
    for c in out:
        out[c].sort()
    return out

def overlap_regions(intervals_a, intervals_b):
    regions = []
    for chrom, a_list in intervals_a.items():
        b_list = intervals_b.get(chrom, [])
        for a_s, a_e in a_list:
            for b_s, b_e in b_list:
                lo, hi = max(a_s, b_s), min(a_e, b_e)
                if lo < hi:
                    regions.append((chrom, lo, hi))
    return regions

def load_raw_calls(bed_path, cnv_type):
    from collections import defaultdict
    out = defaultdict(list)
    if bed_path is None:
        return out
    with open(bed_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            if parts[3] != cnv_type:
                continue
            out[parts[0]].append((int(parts[1]), int(parts[2])))
    return out

def region_overlaps_any(chrom, start, end, intervals_by_chrom):
    for s, e in intervals_by_chrom.get(chrom, []):
        if max(start, s) < min(end, e):
            return True
    return False

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc-log", required=True)
    ap.add_argument("--merge-dir", required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--resolution", default="25000")
    ap.add_argument("--min-support", type=int, default=1,
                     help="qc log min_support level . "
                          "min_support=1( ) conflict .")
    ap.add_argument("--out-dir", default="./conflict_diagnosis")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # QC (population, sample) 
    targets = []
    with open(args.qc_log) as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            pop, sample, min_sup = parts[0], parts[1], int(parts[2])
            if min_sup == args.min_support:
                targets.append((pop, sample))

    print(f" : {len(targets)} (population, sample) "
          f"(min_support={args.min_support})", file=sys.stderr)

    combo_counter = Counter()
    tool_involvement = Counter() # tool -> (DEL , DUP ) 
    tool_del_count = Counter()
    tool_dup_count = Counter()
    detail_path = os.path.join(args.out_dir, "conflict_detail.csv")

    with open(detail_path, "w", newline="") as detail_f:
        writer = csv.writer(detail_f)
        writer.writerow(["population", "sample", "chrom", "start", "end",
                          "del_supporting_tools", "dup_supporting_tools"])

        for pop, sample in targets:
            del_merge = find_merge_bed(args.merge_dir, sample, args.resolution, "DEL")
            dup_merge = find_merge_bed(args.merge_dir, sample, args.resolution, "DUP")
            if del_merge is None or dup_merge is None:
                continue

            del_presence = load_presence_intervals(del_merge, args.min_support)
            dup_presence = load_presence_intervals(dup_merge, args.min_support)
            regions = overlap_regions(del_presence, dup_presence)
            if not regions:
                continue

            # 4 tool raw 
            tool_del_calls = {}
            tool_dup_calls = {}
            for tool in TOOLS:
                raw = find_raw_bed(args.data_dir, sample, tool, args.resolution)
                tool_del_calls[tool] = load_raw_calls(raw, "DEL")
                tool_dup_calls[tool] = load_raw_calls(raw, "DUP")

            for chrom, start, end in regions:
                del_tools = tuple(sorted(
                    t for t in TOOLS if region_overlaps_any(chrom, start, end, tool_del_calls[t])
                ))
                dup_tools = tuple(sorted(
                    t for t in TOOLS if region_overlaps_any(chrom, start, end, tool_dup_calls[t])
                ))
                if not del_tools or not dup_tools:
                    # raw ( ) 
                    continue
                writer.writerow([pop, sample, chrom, start, end,
                                  ",".join(del_tools), ",".join(dup_tools)])
                combo_counter[(del_tools, dup_tools)] += 1
                for t in del_tools:
                    tool_del_count[t] += 1
                for t in dup_tools:
                    tool_dup_count[t] += 1

    combo_path = os.path.join(args.out_dir, "conflict_combo_summary.csv")
    with open(combo_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["del_tools", "dup_tools", "n_occurrences"])
        for (del_tools, dup_tools), n in combo_counter.most_common():
            writer.writerow([",".join(del_tools), ",".join(dup_tools), n])

    tool_path = os.path.join(args.out_dir, "conflict_tool_involvement.csv")
    with open(tool_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["tool", "n_times_on_DEL_side", "n_times_on_DUP_side", "n_times_total"])
        for tool in TOOLS:
            nd, nu = tool_del_count[tool], tool_dup_count[tool]
            writer.writerow([tool, nd, nu, nd + nu])

    print(f"\nSaved to: {detail_path}")
    print(f" : {combo_path}")
    print(f"tool : {tool_path}")
    print("\n=== Tool conflict ( ) ===")
    for tool in sorted(TOOLS, key=lambda t: -(tool_del_count[t] + tool_dup_count[t])):
        print(f" {tool}: DEL {tool_del_count[tool]}, DUP {tool_dup_count[tool]}, "
              f" {tool_del_count[tool]+tool_dup_count[tool]}")

if __name__ == "__main__":
    main()
