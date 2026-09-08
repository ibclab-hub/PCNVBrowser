#!/usr/bin/env python3
# Role: Permutation test (bedtools shuffle) for the empirical significance of population-level concordance (Supplementary Table S8).
import argparse
import csv
import gzip
import os
import random
import re
import subprocess
import sys
import tempfile

SIZE_BINS = [(50, 500), (500, 1000), (1000, 5000), (5000, 25000), (25000, float("inf"))]

def bin_label(lo, hi):
    if hi == float("inf"):
        return f">={lo}bp"
    return f"{lo}-{hi}bp"

def size_bin_of(size):
    for lo, hi in SIZE_BINS:
        if lo <= size < hi:
            return (lo, hi)
    return None

def load_bed(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            chrom, start, end, svtype, freq = parts[:5]
            try:
                start = int(start)
                end = int(end)
                freq = float(freq)
            except ValueError:
                continue
            if end <= start:
                continue
            rows.append((chrom, start, end, svtype, freq))
    return rows

def generate_genome_from_vcf(vcf_path, out_path):
    opener = gzip.open if vcf_path.endswith(".gz") else open
    pattern = re.compile(r"ID=([^,]+),length=(\d+)")
    contigs = []
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                break
            if line.startswith("##contig"):
                m = pattern.search(line)
                if m:
                    contigs.append((m.group(1), int(m.group(2))))
    if not contigs:
        sys.exit(f"{vcf_path} ##contig ")
    with open(out_path, "w") as out:
        for chrom, length in contigs:
            out.write(f"{chrom}\t{length}\n")
    print(f"genome : {out_path} ({len(contigs)} contig)", file=sys.stderr)
    return out_path

def write_our_indexed_bed(rows, path):
    idx_to_size = {}
    idx_to_bin = {}
    with open(path, "w") as f:
        for i, (chrom, start, end, _t, _f) in enumerate(rows):
            f.write(f"{chrom}\t{start}\t{end}\t{i}\n")
            size = end - start
            idx_to_size[i] = size
            idx_to_bin[i] = size_bin_of(size)
    return idx_to_bin

def write_kgp_bed(rows, path):
    with open(path, "w") as f:
        for chrom, start, end, _t, _f in sorted(rows, key=lambda r: (r[0], r[1], r[2])):
            f.write(f"{chrom}\t{start}\t{end}\n")

def matched_idx_set(a_indexed_bed, b_bed, frac):
    try:
        result = subprocess.run(
            ["bedtools", "intersect", "-a", a_indexed_bed, "-b", b_bed,
             "-f", str(frac), "-r", "-u"],
            capture_output=True, text=True, check=True,
        )
    except subprocess.CalledProcessError as e:
        print(f" [] bedtools intersect : {e.stderr}", file=sys.stderr)
        return set()
    idxs = set()
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        if len(f) >= 4:
            idxs.add(int(f[3]))
    return idxs

def rates_by_bin(idx_set, idx_to_bin):
    bin_total = {b: 0 for b in SIZE_BINS}
    bin_matched = {b: 0 for b in SIZE_BINS}
    for idx, b in idx_to_bin.items():
        if b is None:
            continue
        bin_total[b] += 1
        if idx in idx_set:
            bin_matched[b] += 1
    return bin_total, bin_matched

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--our-dir", required=True)
    ap.add_argument("--kgp-dir", required=True)
    ap.add_argument("--resolution", default="5000")
    ap.add_argument("--populations", required=True)
    ap.add_argument("--types", default="DUP,DEL")
    ap.add_argument("--our-pattern", default="{pop}.tools_{res}.{type}.carrier_freq_min1.bed")
    ap.add_argument("--kgp-pattern", default="{pop}.1kgp_sv.{type}.freq.bed")
    ap.add_argument("--min-size", type=int, default=0,
                     help="(bp) (compare_1kgp_concordance.py )")
    ap.add_argument("--reciprocal-frac", type=float, default=0.5)
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--genome", default=None, help="chrom\\tsize ( )")
    ap.add_argument("--vcf", default=None,
                     help="--genome , VCF ##contig genome ")
    ap.add_argument("--keep-chrom", action="store_true",
                     help="shuffle (bedtools shuffle -chrom). "
                          " genome shuffle")
    ap.add_argument("--exclude", default=None,
                     help="shuffle (: assembly gap/centromere) bed . ")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out-dir", default="./tier2_permutation_results")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    if args.genome:
        genome_path = args.genome
        if not os.path.isfile(genome_path):
            sys.exit(f"genome : {genome_path}")
    elif args.vcf:
        genome_path = os.path.join(args.out_dir, "genome.txt")
        generate_genome_from_vcf(args.vcf, genome_path)
    else:
        sys.exit("--genome --vcf ")

    populations = [p.strip() for p in args.populations.split(",") if p.strip()]
    types = [t.strip() for t in args.types.split(",") if t.strip()]

    out_path = os.path.join(args.out_dir, "permutation_null_by_size.csv")
    missing_log_path = os.path.join(args.out_dir, "missing_files.log")

    random.seed(args.seed)

    with open(out_path, "w", newline="") as fo, \
         open(missing_log_path, "w") as fm, \
         tempfile.TemporaryDirectory() as tmpdir:

        writer = csv.writer(fo)
        writer.writerow(["population", "type", "size_bin", "n_in_bin",
                          "observed_rate", "null_mean", "null_sd",
                          "z_score", "empirical_p", "n_perm"])

        for pop in populations:
            for svtype in types:
                our_path = os.path.join(
                    args.our_dir, args.our_pattern.format(pop=pop, res=args.resolution, type=svtype))
                kgp_path = os.path.join(
                    args.kgp_dir, args.kgp_pattern.format(pop=pop, type=svtype))

                if not os.path.isfile(our_path) or not os.path.isfile(kgp_path):
                    fm.write(f"{pop}\t{svtype}\tFILE_MISSING\tour={our_path} kgp={kgp_path}\n")
                    continue

                our_rows_raw = load_bed(our_path)
                kgp_rows = load_bed(kgp_path)
                if args.min_size > 0:
                    our_rows = [r for r in our_rows_raw if (r[2] - r[1]) >= args.min_size]
                else:
                    our_rows = our_rows_raw
                if not our_rows or not kgp_rows:
                    fm.write(f"{pop}\t{svtype}\tEMPTY_AFTER_FILTER\n")
                    continue

                our_indexed_path = os.path.join(tmpdir, f"{pop}_{svtype}_our_idx.bed")
                kgp_bed_path = os.path.join(tmpdir, f"{pop}_{svtype}_kgp.bed")
                idx_to_bin = write_our_indexed_bed(our_rows, our_indexed_path)
                write_kgp_bed(kgp_rows, kgp_bed_path)

                observed_idx = matched_idx_set(our_indexed_path, kgp_bed_path, args.reciprocal_frac)
                obs_total, obs_matched = rates_by_bin(observed_idx, idx_to_bin)

                # null distribution: n_perm shuffle
                null_rates = {b: [] for b in SIZE_BINS}
                shuffle_cmd_base = ["bedtools", "shuffle", "-i", our_indexed_path,
                                     "-g", genome_path]
                if args.keep_chrom:
                    shuffle_cmd_base.append("-chrom")
                if args.exclude:
                    shuffle_cmd_base += ["-excl", args.exclude]

                for p in range(args.n_perm):
                    shuffled_path = os.path.join(tmpdir, f"{pop}_{svtype}_shuf_{p}.bed")
                    seed_i = args.seed + p + 1
                    cmd = shuffle_cmd_base + ["-seed", str(seed_i)]
                    try:
                        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                    except subprocess.CalledProcessError as e:
                        print(f" [] {pop} {svtype} perm {p} shuffle : {e.stderr}",
                              file=sys.stderr)
                        continue
                    with open(shuffled_path, "w") as sf:
                        sf.write(result.stdout)

                    shuf_idx = matched_idx_set(shuffled_path, kgp_bed_path, args.reciprocal_frac)
                    _, shuf_matched = rates_by_bin(shuf_idx, idx_to_bin)
                    for b in SIZE_BINS:
                        if obs_total[b] > 0:
                            null_rates[b].append(shuf_matched[b] / obs_total[b])

                for b in SIZE_BINS:
                    label = bin_label(*b)
                    n_total = obs_total[b]
                    if n_total == 0:
                        writer.writerow([pop, svtype, label, 0, "NA", "NA", "NA",
                                          "NA", "NA", args.n_perm])
                        continue
                    observed_rate = obs_matched[b] / n_total
                    nulls = null_rates[b]
                    if not nulls:
                        writer.writerow([pop, svtype, label, n_total,
                                          f"{observed_rate:.4f}", "NA", "NA", "NA", "NA",
                                          args.n_perm])
                        continue
                    null_mean = sum(nulls) / len(nulls)
                    null_var = sum((x - null_mean) ** 2 for x in nulls) / max(len(nulls) - 1, 1)
                    null_sd = null_var ** 0.5
                    z = (observed_rate - null_mean) / null_sd if null_sd > 0 else float("nan")
                    n_geq = sum(1 for x in nulls if x >= observed_rate)
                    emp_p = (n_geq + 1) / (len(nulls) + 1)

                    writer.writerow([
                        pop, svtype, label, n_total,
                        f"{observed_rate:.4f}", f"{null_mean:.4f}", f"{null_sd:.4f}",
                        f"{z:.2f}" if z == z else "NA",
                        f"{emp_p:.4f}", args.n_perm,
                    ])

                print(f"{pop} {svtype} ({args.n_perm} permutation)", file=sys.stderr)

    print(f"\n. : {out_path}")
    print(f" : {missing_log_path}")

if __name__ == "__main__":
    main()
