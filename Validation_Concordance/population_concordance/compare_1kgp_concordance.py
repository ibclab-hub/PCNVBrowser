#!/usr/bin/env python3
# Role: Compares PCNVBrowser population-level frequencies (carrier frequency and call-set support) with the 1000 Genomes SV callset (Supplementary Table S8).
import argparse
import csv
import os
import subprocess
import sys
import tempfile

try:
    from scipy import stats as _stats
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

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

def write_bed_sorted(rows, path):
    rows_sorted = sorted(rows, key=lambda r: (r[0], r[1], r[2]))
    with open(path, "w") as f:
        for chrom, start, end, svtype, freq in rows_sorted:
            f.write(f"{chrom}\t{start}\t{end}\t{svtype}\t{freq}\n")
    return len(rows_sorted)

def reciprocal_match(our_path, kgp_path, frac):
    try:
        result = subprocess.run(
            ["bedtools", "intersect", "-a", our_path, "-b", kgp_path,
             "-f", str(frac), "-r", "-wo"],
            capture_output=True, text=True, check=True,
        )
    except subprocess.CalledProcessError as e:
        print(f" [] bedtools intersect : {e.stderr}", file=sys.stderr)
        return {}
    except FileNotFoundError:
        sys.exit("bedtools . PATH .")

    best = {}
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        # our: chrom,start,end,type,freq (5) / kgp: chrom,start,end,type,freq (5) / overlap_bp (1)
        if len(f) < 11:
            continue
        our_chrom, our_start, our_end = f[0], int(f[1]), int(f[2])
        our_freq = float(f[4])
        kgp_chrom, kgp_start, kgp_end = f[5], int(f[6]), int(f[7])
        kgp_freq = float(f[9])
        overlap_bp = int(f[10])

        key = (our_chrom, our_start, our_end)
        if key not in best or overlap_bp > best[key][4]:
            best[key] = (kgp_chrom, kgp_start, kgp_end, kgp_freq, overlap_bp)
    return best

def correlation(xs, ys):
    if len(xs) < 3:
        return None, None, None, None
    if HAVE_SCIPY:
        pr, pp = _stats.pearsonr(xs, ys)
        sr, sp = _stats.spearmanr(xs, ys)
        return pr, pp, sr, sp
    return None, None, None, None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--our-dir", required=True)
    ap.add_argument("--kgp-dir", required=True)
    ap.add_argument("--resolution", default="5000",
                     help="resolution (: 5000)")
    ap.add_argument("--populations", required=True,
                     help="population ( )")
    ap.add_argument("--types", default="DUP,DEL",
                     help="SV (: DUP,DEL)")
    ap.add_argument("--our-pattern", default="{pop}.tools_{res}.{type}.carrier_freq_min1.bed",
                     help="({pop},{res},{type} )")
    ap.add_argument("--kgp-pattern", default="{pop}.1kgp_sv.{type}.freq.bed",
                     help="1000G ({pop},{type} )")
    ap.add_argument("--reciprocal-frac", type=float, default=0.5)
    ap.add_argument("--min-size", type=int, default=0,
                     help="region (bp) . "
                          "PCNVBrowser calling resolution(: 5000) , "
                          "≥1bp overlap merge sub-resolution (sliver) "
                          "(0= , )")
    ap.add_argument("--out-dir", default="./tier2_concordance_results")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    populations = [p.strip() for p in args.populations.split(",") if p.strip()]
    types = [t.strip() for t in args.types.split(",") if t.strip()]

    concordance_path = os.path.join(args.out_dir, "concordance_by_size.csv")
    correlation_path = os.path.join(args.out_dir, "af_correlation.csv")
    missing_log_path = os.path.join(args.out_dir, "missing_files.log")
    filtered_log_path = os.path.join(args.out_dir, "min_size_filtered.csv")

    with open(concordance_path, "w", newline="") as fc, \
         open(correlation_path, "w", newline="") as fr, \
         open(missing_log_path, "w") as fm, \
         open(filtered_log_path, "w", newline="") as ff, \
         tempfile.TemporaryDirectory() as tmpdir:

        writer_f = csv.writer(ff)
        writer_f.writerow(["population", "type", "min_size_bp",
                            "n_our_total", "n_our_filtered_out", "n_our_kept"])

        writer_c = csv.writer(fc)
        writer_c.writerow(["population", "type", "size_bin", "n_our_in_bin",
                            "n_matched", "concordance_rate"])
        writer_r = csv.writer(fr)
        writer_r.writerow(["population", "type", "n_pairs",
                            "pearson_r", "pearson_p", "spearman_r", "spearman_p"])

        for pop in populations:
            for svtype in types:
                our_path = os.path.join(
                    args.our_dir,
                    args.our_pattern.format(pop=pop, res=args.resolution, type=svtype))
                kgp_path = os.path.join(
                    args.kgp_dir,
                    args.kgp_pattern.format(pop=pop, type=svtype))

                if not os.path.isfile(our_path):
                    fm.write(f"{pop}\t{svtype}\tOUR_FILE_MISSING\t{our_path}\n")
                    continue
                if not os.path.isfile(kgp_path):
                    fm.write(f"{pop}\t{svtype}\tKGP_FILE_MISSING\t{kgp_path}\n")
                    continue

                our_rows_raw = load_bed(our_path)
                kgp_rows = load_bed(kgp_path)

                # min-size : PCNVBrowser calling resolution region
                # ≥1bp overlap merge sub-resolution (sliver) 
                # . sensitivity
                # analysis .
                if args.min_size > 0:
                    n_before = len(our_rows_raw)
                    our_rows = [r for r in our_rows_raw if (r[2] - r[1]) >= args.min_size]
                    n_filtered = n_before - len(our_rows)
                    writer_f.writerow([pop, svtype, args.min_size,
                                        n_before, n_filtered, len(our_rows)])
                else:
                    our_rows = our_rows_raw

                if not our_rows or not kgp_rows:
                    fm.write(f"{pop}\t{svtype}\tEMPTY_AFTER_PARSE\t"
                             f"our={len(our_rows)} kgp={len(kgp_rows)}\n")
                    continue

                our_sorted_path = os.path.join(tmpdir, f"{pop}_{svtype}_our.bed")
                kgp_sorted_path = os.path.join(tmpdir, f"{pop}_{svtype}_kgp.bed")
                write_bed_sorted(our_rows, our_sorted_path)
                write_bed_sorted(kgp_rows, kgp_sorted_path)

                matches = reciprocal_match(our_sorted_path, kgp_sorted_path,
                                            args.reciprocal_frac)

                # (1) concordance rate: region 
                bin_totals = {b: 0 for b in SIZE_BINS}
                bin_matched = {b: 0 for b in SIZE_BINS}
                for chrom, start, end, _t, _f in our_rows:
                    b = size_bin_of(end - start)
                    if b is None:
                        continue
                    bin_totals[b] += 1
                    if (chrom, start, end) in matches:
                        bin_matched[b] += 1

                for b in SIZE_BINS:
                    label = bin_label(*b)
                    n_total = bin_totals[b]
                    n_matched = bin_matched[b]
                    rate = f"{n_matched / n_total:.4f}" if n_total > 0 else "NA"
                    writer_c.writerow([pop, svtype, label, n_total, n_matched, rate])

                # (2) AF correlation
                our_freq_by_key = {(c, s, e): fr_ for c, s, e, _t, fr_ in our_rows}
                xs, ys = [], []
                for key, (_kc, _ks, _ke, kgp_freq, _bp) in matches.items():
                    xs.append(our_freq_by_key[key])
                    ys.append(kgp_freq)

                pr, pp, sr, sp = correlation(xs, ys)
                writer_r.writerow([
                    pop, svtype, len(xs),
                    f"{pr:.4f}" if pr is not None else "NA",
                    f"{pp:.2e}" if pp is not None else "NA",
                    f"{sr:.4f}" if sr is not None else "NA",
                    f"{sp:.2e}" if sp is not None else "NA",
                ])

                print(f"{pop} {svtype}: our={len(our_rows)} kgp={len(kgp_rows)} "
                      f"matched={len(matches)}", file=sys.stderr)

    if not HAVE_SCIPY:
        print("\n[] scipy — correlation NA . "
              "pip install scipy .", file=sys.stderr)

    print(f"\n.")
    print(f"concordance (): {concordance_path}")
    print(f"AF correlation: {correlation_path}")
    print(f" : {missing_log_path}")

if __name__ == "__main__":
    main()
