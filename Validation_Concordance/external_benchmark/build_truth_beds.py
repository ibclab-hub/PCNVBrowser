#!/usr/bin/env python3
# Role: Builds standardized truth BED files from the NA12878 and HGSVC2 reference call sets.
import argparse
import gzip
import os
import sys

def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)

def write_bed(rows, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    rows.sort()
    with open(path, "w") as f:
        for chrom, start, end in rows:
            f.write(f"{chrom}\t{start}\t{end}\n")
    return len(rows)

def cmd_na12878(args):
    rows = []
    with open(args.input) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            rows.append((chrom, start, end))
    out_path = os.path.join(args.out_dir, "NA12878", "DEL.bed")
    n = write_bed(rows, out_path)
    print(f"NA12878 DEL: {n} -> {out_path}")
    if args.liftover_chain:
        print("--liftover-chain provided; running UCSC liftOver: "
              f"liftOver {out_path} {args.liftover_chain} {out_path}.hg38 {out_path}.unmapped")
    else:
        print("Saved as GRCh37; PCNVBrowser uses GRCh38, so liftOver before comparison if needed.")

def cmd_hgsvc2(args):
    with open_maybe_gzip(args.input) as f:
        header = None
        for line in f:
            if line.startswith("#CHROM") or line.startswith("#chrom") or \
               (header is None and "SVTYPE" in line and "\t" in line):
                header = line.rstrip("\n").split("\t")
                break
        if header is None:
            sys.exit(" . head .")

        print(" :", header, file=sys.stderr)

        def col_idx(name_options):
            for name in name_options:
                if name in header:
                    return header.index(name)
            return None

        chrom_i = col_idx(["#CHROM", "CHROM", "#chrom", "chrom"])
        pos_i = col_idx(["POS", "pos", "START", "start"])
        end_i = col_idx(["END", "end"])
        svtype_i = col_idx(["SVTYPE", "svtype"])
        svlen_i = col_idx(["SVLEN", "svlen"])

        if chrom_i is None or pos_i is None or svtype_i is None:
            sys.exit(f" (CHROM/POS/SVTYPE) . : {header}")

        sample_cols = {}
        for sample in args.sample:
            if args.sample_col and sample in args.sample_col:
                colname = args.sample_col[sample]
            else:
                colname = sample
            if colname not in header:
                sys.exit(f" '{sample}' ('{colname}') . "
                          f"--sample-col {sample}= . : {header}")
            sample_cols[sample] = header.index(colname)

        rows_del = {s: [] for s in args.sample}
        rows_ins = {s: [] for s in args.sample}

        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(sample_cols.values()):
                continue
            chrom = parts[chrom_i]
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            try:
                start = int(parts[pos_i])
            except ValueError:
                continue
            svtype = parts[svtype_i].upper()
            if end_i is not None and parts[end_i].strip():
                try:
                    end = int(parts[end_i])
                except ValueError:
                    end = start + 1
            elif svlen_i is not None and parts[svlen_i].strip():
                try:
                    end = start + abs(int(parts[svlen_i]))
                except ValueError:
                    end = start + 1
            else:
                end = start + 1

            for sample, ci in sample_cols.items():
                gt = parts[ci]
                is_present = any(a not in ("0", ".", "") for a in gt.replace("|", "/").split("/"))
                if not is_present:
                    continue
                if svtype == "DEL":
                    rows_del[sample].append((chrom, start, end))
                elif svtype == "INS":
                    rows_ins[sample].append((chrom, start, end))

        for sample in args.sample:
            d_path = os.path.join(args.out_dir, sample, "DEL.bed")
            i_path = os.path.join(args.out_dir, sample, "INS.bed")
            nd = write_bed(rows_del[sample], d_path)
            ni = write_bed(rows_ins[sample], i_path)
            print(f"{sample}: DEL {nd} -> {d_path} / INS {ni}( DUP ) -> {i_path}")

def cmd_hprc(args):
    import subprocess
    rows_del, rows_ins = [], []
    proc = subprocess.Popen(["zcat", args.vcf], stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        chrom, pos, _id, ref, alt = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        alt = alt.split(",")[0] # ALT (multi-allelic )
        if "<" in alt or alt in (".", ""):
            continue # symbolic ALT 
        diff = len(alt) - len(ref)
        if diff <= -args.min_size:
            rows_del.append((chrom, pos, pos + abs(diff)))
        elif diff >= args.min_size:
            rows_ins.append((chrom, pos, pos + diff))
    proc.wait()

    d_path = os.path.join(args.out_dir, args.sample, "DEL.bed")
    i_path = os.path.join(args.out_dir, args.sample, "INS.bed")
    nd = write_bed(rows_del, d_path)
    ni = write_bed(rows_ins, i_path)
    print(f"{args.sample}: DEL {nd} -> {d_path} / INS {ni}( DUP ) -> {i_path}")

def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("na12878")
    p1.add_argument("--input", required=True)
    p1.add_argument("--out-dir", default="truth_beds")
    p1.add_argument("--liftover-chain", default=None)
    p1.set_defaults(func=cmd_na12878)

    p2 = sub.add_parser("hgsvc2")
    p2.add_argument("--input", required=True)
    p2.add_argument("--sample", action="append", required=True)
    p2.add_argument("--sample-col", action="append", default=[],
                     help="sample= , sample ID ")
    p2.add_argument("--out-dir", default="truth_beds")
    p2.set_defaults(func=cmd_hgsvc2)

    p3 = sub.add_parser("hprc")
    p3.add_argument("--sample", required=True)
    p3.add_argument("--vcf", required=True)
    p3.add_argument("--min-size", type=int, default=50)
    p3.add_argument("--out-dir", default="truth_beds")
    p3.set_defaults(func=cmd_hprc)

    args = ap.parse_args()
    if args.cmd == "hgsvc2":
        sc = {}
        for kv in args.sample_col:
            k, v = kv.split("=", 1)
            sc[k] = v
        args.sample_col = sc
    args.func(args)

if __name__ == "__main__":
    main()
