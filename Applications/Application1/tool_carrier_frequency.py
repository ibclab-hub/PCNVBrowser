#!/usr/bin/env python3
# Role: Computes single-tool (union) carrier frequency per population, used for the caller-specific 'ALL' metrics in Application 1.
import argparse
import json
import os
import subprocess
import sys
import tempfile

TOOLS = ["CN.MOPS", "CNVKIT", "FREEC", "READDEPTH"]
TYPES = ["DEL", "DUP"]

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

def to_type_presence_bed(src_path, cnv_type, dst_path):
    rows = []
    with open(src_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            if parts[3] != cnv_type:
                continue
            rows.append((parts[0], int(parts[1]), int(parts[2])))
    if not rows:
        return None
    rows.sort()
    with open(dst_path, "w") as out:
        for chrom, start, end in rows:
            out.write(f"{chrom}\t{start}\t{end}\n")
    return dst_path

def union_bed(paths, dst_path):
    rows = []
    for p in paths:
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

def write_carrier_freq(presence_beds, n_denom, out_path):
    if not presence_beds:
        return False
    if len(presence_beds) == 1:
        freq = 1 / n_denom
        with open(presence_beds[0]) as f_in, open(out_path, "w") as out_f:
            for line in f_in:
                parts = line.split()
                if len(parts) < 3:
                    continue
                out_f.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\tALL\t{freq:.5f}\n")
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
            out_f.write(f"{chrom}\t{start}\t{end}\tALL\t{freq:.5f}\n")
    return True

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--resolution", default="25000")
    ap.add_argument("--populations", default=None,
                     help="population . ")
    ap.add_argument("--out-dir", default="./tool_carrier_freq_results")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    pop2samples = load_populations(args.metadata)
    target_pops = (args.populations.split(",") if args.populations
                    else list(pop2samples.keys()))

    missing_log_path = os.path.join(args.out_dir, f"missing_files_res{args.resolution}.log")

    with tempfile.TemporaryDirectory() as tmpdir, open(missing_log_path, "w") as miss_f:
        for pop in target_pops:
            samples = pop2samples.get(pop, [])
            if not samples:
                print(f"population '{pop}' , ", file=sys.stderr)
                continue
            n_denom = len(samples)

            for tool in TOOLS:
                presence_beds = []
                for sample in samples:
                    raw = find_raw_bed(args.data_dir, sample, tool, args.resolution)
                    if raw is None:
                        miss_f.write(f"{pop}\t{sample}\t{tool}\tRAW_BED_NOT_FOUND\n")
                        continue
                    del_bed = to_type_presence_bed(
                        raw, "DEL", os.path.join(tmpdir, f"{sample}_{tool}_DEL.bed"))
                    dup_bed = to_type_presence_bed(
                        raw, "DUP", os.path.join(tmpdir, f"{sample}_{tool}_DUP.bed"))
                    union_path = union_bed(
                        [del_bed, dup_bed], os.path.join(tmpdir, f"{sample}_{tool}_ALL.bed"))
                    if union_path:
                        presence_beds.append(union_path)

                out_path = os.path.join(
                    args.out_dir, f"{pop}.{tool}.{args.resolution}.ALL.tool_carrier_freq.bed"
                )
                ok = write_carrier_freq(presence_beds, n_denom, out_path)
                status = f" (carrier bed {len(presence_beds)}, ={n_denom})" if ok else ""
                print(f"{pop} {tool}: {status} -> {out_path}", file=sys.stderr)

    print(f"\n. : {missing_log_path}")

if __name__ == "__main__":
    main()
