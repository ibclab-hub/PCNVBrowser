#!/usr/bin/env python3
# Role: Repeats the 1000 Genomes concordance comparison across reciprocal-overlap thresholds 10-90% (Supplementary Table S9).
import argparse
import csv
import os
import subprocess
import sys

def run_compare(compare_script, common_args, frac, out_subdir):
    cmd = [sys.executable, compare_script] + common_args + [
        "--reciprocal-frac", str(frac),
        "--out-dir", out_subdir,
    ]
    print(f" : frac={frac} -> {out_subdir}", file=sys.stderr)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f" [] frac={frac} :\n{result.stderr}", file=sys.stderr)
        return False
    return True

def concat_csv(paths_with_frac, out_path, extra_col_name="reciprocal_frac"):
    header_written = False
    with open(out_path, "w", newline="") as out_f:
        writer = None
        for frac, path in paths_with_frac:
            if not os.path.isfile(path):
                print(f" [] : {path}", file=sys.stderr)
                continue
            with open(path) as in_f:
                reader = csv.reader(in_f)
                header = next(reader, None)
                if header is None:
                    continue
                if not header_written:
                    writer = csv.writer(out_f)
                    writer.writerow([extra_col_name] + header)
                    header_written = True
                for row in reader:
                    writer.writerow([frac] + row)
    return out_path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--compare-script", required=True,
                     help="compare_1kgp_concordance.py ")
    ap.add_argument("--our-dir", required=True)
    ap.add_argument("--kgp-dir", required=True)
    ap.add_argument("--resolution", default="5000")
    ap.add_argument("--populations", required=True)
    ap.add_argument("--types", default="DUP,DEL")
    ap.add_argument("--our-pattern", default=None,
                     help="compare_1kgp_concordance.py --our-pattern ()")
    ap.add_argument("--kgp-pattern", default=None,
                     help="compare_1kgp_concordance.py --kgp-pattern ()")
    ap.add_argument("--min-size", type=int, default=0)
    ap.add_argument("--fracs", default="0.1,0.3,0.5,0.7,0.9",
                     help="reciprocal-overlap fraction ")
    ap.add_argument("--out-dir", default="./tier2_frac_sweep")
    args = ap.parse_args()

    if not os.path.isfile(args.compare_script):
        sys.exit(f"compare-script : {args.compare_script}")

    os.makedirs(args.out_dir, exist_ok=True)
    fracs = [f.strip() for f in args.fracs.split(",") if f.strip()]

    common_args = [
        "--our-dir", args.our_dir,
        "--kgp-dir", args.kgp_dir,
        "--resolution", args.resolution,
        "--populations", args.populations,
        "--types", args.types,
        "--min-size", str(args.min_size),
    ]
    if args.our_pattern:
        common_args += ["--our-pattern", args.our_pattern]
    if args.kgp_pattern:
        common_args += ["--kgp-pattern", args.kgp_pattern]

    concordance_pairs = []
    correlation_pairs = []

    for frac in fracs:
        out_subdir = os.path.join(args.out_dir, f"frac_{frac}")
        os.makedirs(out_subdir, exist_ok=True)
        ok = run_compare(args.compare_script, common_args, frac, out_subdir)
        if not ok:
            continue
        concordance_pairs.append((frac, os.path.join(out_subdir, "concordance_by_size.csv")))
        correlation_pairs.append((frac, os.path.join(out_subdir, "af_correlation.csv")))

    concordance_out = os.path.join(args.out_dir, "concordance_sweep.csv")
    correlation_out = os.path.join(args.out_dir, "af_correlation_sweep.csv")
    concat_csv(concordance_pairs, concordance_out)
    concat_csv(correlation_pairs, correlation_out)

    print(f"\n.")
    print(f"concordance sweep (frac): {concordance_out}")
    print(f"AF correlation sweep (frac): {correlation_out}")
    print(f" frac : {args.out_dir}/frac_<value>/")

if __name__ == "__main__":
    main()
