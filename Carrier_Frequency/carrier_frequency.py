#!/usr/bin/env python3
# Role: Computes individual-based carrier frequency (DEL/DUP/ALL, support thresholds >=1-4) per population from integrated CNV calls.
import argparse
import json
import os
import subprocess
import sys
import tempfile

TYPES = ["DEL", "DUP"]
MIN_SUPPORT_LEVELS = (1, 2, 3, 4)
FRAC_TO_NTOOLS = {0.25: 1, 0.5: 2, 0.75: 3, 1.0: 4}

def load_populations(metadata_path):
    with open(metadata_path) as f:
        meta = json.load(f)
    pop2samples = {}
    for pop, samples in meta.items():
        pop2samples[pop] = [s for s in samples if not s.startswith("merged_")]
    return pop2samples

def find_merge_bed(merge_dir, sample, resolution, cnv_type):
    path = os.path.join(merge_dir, f"{sample}.tools_{resolution}.{cnv_type}.merge_sort.bed")
    return path if os.path.isfile(path) else None

def to_presence_bed(src_path, dst_path, min_support):
    rows = []
    with open(src_path) as f:
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
            rows.append((parts[0], int(parts[1]), int(parts[2])))
    if not rows:
        return None
    rows.sort()
    with open(dst_path, "w") as out:
        for chrom, start, end in rows:
            out.write(f"{chrom}\t{start}\t{end}\n")
    return dst_path

def count_overlap(bed_a, bed_b):
    try:
        result = subprocess.run(
            ["bedtools", "intersect", "-a", bed_a, "-b", bed_b],
            capture_output=True, text=True, check=True,
        )
    except subprocess.CalledProcessError:
        return 0, 0
    n_intervals = 0
    bp = 0
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        n_intervals += 1
        bp += int(parts[2]) - int(parts[1])
    return n_intervals, bp

def union_bed(bed_a, bed_b, dst_path):
    rows = []
    for p in (bed_a, bed_b):
        if p is None:
            continue
        with open(p) as f:
            for line in f:
                if line.strip():
                    rows.append(line.strip())
    if not rows:
        return None
    cat_path = dst_path + ".cat"
    with open(cat_path, "w") as f:
        for r in sorted(set(rows)):
            f.write(r + "\n")
    try:
        subprocess.run(f"sort -k1,1 -k2,2n {cat_path} | bedtools merge -i - > {dst_path}",
                        shell=True, check=True)
    except subprocess.CalledProcessError:
        return None
    os.remove(cat_path)
    return dst_path if os.path.getsize(dst_path) > 0 else None

def write_carrier_freq(presence_beds, n_denom, cnv_type_label, out_path):
    if not presence_beds:
        return False
    if len(presence_beds) == 1:
        freq = 1 / n_denom
        with open(presence_beds[0]) as f_in, open(out_path, "w") as out_f:
            for line in f_in:
                parts = line.split()
                if len(parts) < 3:
                    continue
                out_f.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{cnv_type_label}\t{freq:.5f}\n")
        return True

    try:
        result = subprocess.run(
            ["bedtools", "multiinter", "-i"] + presence_beds,
            capture_output=True, text=True, check=True,
        )
    except subprocess.CalledProcessError:
        return False

    with open(out_path, "w") as out_f:
        for line in result.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end, n_carrier = parts[0], parts[1], parts[2], parts[3]
            freq = int(n_carrier) / n_denom
            out_f.write(f"{chrom}\t{start}\t{end}\t{cnv_type_label}\t{freq:.5f}\n")
    return True

def load_intervals(path):
    from collections import defaultdict
    out = defaultdict(list)
    if not path or not os.path.isfile(path):
        return out
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            out[parts[0]].append((int(parts[1]), int(parts[2])))
    for c in out:
        out[c].sort()
    return out

def overlaps_any(chrom, start, end, intervals_by_chrom):
    for s, e in intervals_by_chrom.get(chrom, []):
        if max(start, s) < min(end, e):
            return True
    return False

def relabel_all_by_source(all_path, del_path, dup_path):
    if not os.path.isfile(all_path):
        return
    del_intervals = load_intervals(del_path)
    dup_intervals = load_intervals(dup_path)

    rows = []
    with open(all_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 5:
                continue
            chrom, start, end, _old_label, freq = parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4]
            has_del = overlaps_any(chrom, start, end, del_intervals)
            has_dup = overlaps_any(chrom, start, end, dup_intervals)
            if has_del and has_dup:
                label = "ALL"
            elif has_del:
                label = "DEL"
            elif has_dup:
                label = "DUP"
            else:
                label = "ALL" # DEL/DUP ALL -> 
            rows.append((chrom, start, end, label, freq))

    with open(all_path, "w") as f:
        for chrom, start, end, label, freq in rows:
            f.write(f"{chrom}\t{start}\t{end}\t{label}\t{freq}\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--merge-dir", required=True,
                     help="{sample}.tools_{res}.{DEL|DUP}.merge_sort.bed ")
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--resolution", default="25000")
    ap.add_argument("--out-dir", default="./carrier_freq_results")
    ap.add_argument("--populations", default=None,
                     help="population . metadata population ")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    pop2samples = load_populations(args.metadata)
    target_pops = (args.populations.split(",") if args.populations
                    else list(pop2samples.keys()))

    missing_log_path = os.path.join(args.out_dir, f"missing_files_res{args.resolution}.log")
    overlap_log_path = os.path.join(args.out_dir, f"del_dup_overlap_qc_res{args.resolution}.log")

    with tempfile.TemporaryDirectory() as tmpdir, \
         open(missing_log_path, "w") as miss_f, \
         open(overlap_log_path, "w") as ov_f:

        ov_f.write("population\tsample\tmin_support\tn_overlap_intervals\toverlap_bp\n")

        for pop in target_pops:
            samples = pop2samples.get(pop, [])
            if not samples:
                print(f"population '{pop}' , ", file=sys.stderr)
                continue
            n_denom = len(samples)

            for min_support in MIN_SUPPORT_LEVELS:
                # 1) DEL / DUP presence bed ()
                del_beds, dup_beds, all_beds = [], [], []

                for sample in samples:
                    src_del = find_merge_bed(args.merge_dir, sample, args.resolution, "DEL")
                    src_dup = find_merge_bed(args.merge_dir, sample, args.resolution, "DUP")
                    if src_del is None and min_support == MIN_SUPPORT_LEVELS[0]:
                        miss_f.write(f"{pop}\t{sample}\tDEL\tMERGE_FILE_NOT_FOUND\n")
                    if src_dup is None and min_support == MIN_SUPPORT_LEVELS[0]:
                        miss_f.write(f"{pop}\t{sample}\tDUP\tMERGE_FILE_NOT_FOUND\n")

                    pres_del = pres_dup = None
                    if src_del is not None:
                        dst = os.path.join(tmpdir, f"{sample}_DEL_min{min_support}.bed")
                        pres_del = to_presence_bed(src_del, dst, min_support)
                        if pres_del:
                            del_beds.append(pres_del)
                    if src_dup is not None:
                        dst = os.path.join(tmpdir, f"{sample}_DUP_min{min_support}.bed")
                        pres_dup = to_presence_bed(src_dup, dst, min_support)
                        if pres_dup:
                            dup_beds.append(pres_dup)

                    # QC: DEL-DUP 
                    if pres_del and pres_dup:
                        n_ov, bp_ov = count_overlap(pres_del, pres_dup)
                        if n_ov > 0:
                            ov_f.write(f"{pop}\t{sample}\t{min_support}\t{n_ov}\t{bp_ov}\n")

                    # ALL(any type) = DEL presence ∪ DUP presence
                    if pres_del or pres_dup:
                        dst_all = os.path.join(tmpdir, f"{sample}_ALL_min{min_support}.bed")
                        pres_all = union_bed(pres_del, pres_dup, dst_all)
                        if pres_all:
                            all_beds.append(pres_all)

                # 2) DEL / DUP / ALL population carrier frequency 
                out_paths = {}
                for type_label, beds in (("DEL", del_beds), ("DUP", dup_beds), ("ALL", all_beds)):
                    out_path = os.path.join(
                        args.out_dir,
                        f"{pop}.tools_{args.resolution}.{type_label}.carrier_freq_min{min_support}.bed"
                    )
                    ok = write_carrier_freq(beds, n_denom, type_label, out_path)
                    out_paths[type_label] = out_path if ok else None
                    status = f" (carrier bed {len(beds)}, ={n_denom})" if ok else "( )"
                    print(f"{pop} {type_label} min{min_support}: {status} -> {out_path}",
                          file=sys.stderr)

                # 3) ALL 4 DEL/DUP 
                # (DEL/DUP/MIXED) ( 5 )
                if out_paths.get("ALL"):
                    relabel_all_by_source(out_paths["ALL"], out_paths.get("DEL"), out_paths.get("DUP"))
                    print(f"{pop} ALL min{min_support}: DEL/DUP/MIXED ",
                          file=sys.stderr)

    print(f"\n. : {missing_log_path}")
    print(f"DEL-DUP QC : {overlap_log_path} ( 0)")

if __name__ == "__main__":
    main()
