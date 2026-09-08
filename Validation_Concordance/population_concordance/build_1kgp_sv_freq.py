#!/usr/bin/env python3
# Role: Extracts population-level DEL/DUP carrier frequency from the 1000 Genomes high-coverage SV callset.
import argparse
import gzip
import json
import os
import subprocess
import sys

def load_populations(metadata_path, wanted=None):
    with open(metadata_path) as f:
        meta = json.load(f)
    pop2samples = {}
    for pop, samples in meta.items():
        if wanted and pop not in wanted:
            continue
        pop2samples[pop] = [s for s in samples if not s.startswith("merged_")]
    return pop2samples

def parse_info(info_str):
    out = {}
    for token in info_str.split(";"):
        if "=" in token:
            k, v = token.split("=", 1)
            out[k] = v
        else:
            out[token] = True
    return out

def open_vcf(path):
    if path.endswith(".gz"):
        # zcat gzip (3202- )
        proc = subprocess.Popen(["zcat", path], stdout=subprocess.PIPE, text=True,
                                 bufsize=1 << 20)
        return proc.stdout, proc
    return open(path), None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--out-dir", default="./tier2_sv_freq")
    ap.add_argument("--populations", default=None,
                     help="population . metadata ")
    ap.add_argument("--include-del", action="store_true",
                     help="DUP DEL ")
    ap.add_argument("--require-pass", action="store_true", default=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    wanted = set(args.populations.split(",")) if args.populations else None
    pop2samples = load_populations(args.metadata, wanted)
    print(f" population {len(pop2samples)}: {list(pop2samples.keys())}", file=sys.stderr)

    stream, proc = open_vcf(args.vcf)

    header = None
    sample_col_idx = {}  # sample_name -> column index (0-based, in full line.split("\t"))
    for line in stream:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            header = line.rstrip("\n").split("\t")
            for i, col in enumerate(header):
                if i >= 9:
                    sample_col_idx[col] = i
            break
    if header is None:
        sys.exit("VCF (#CHROM ) ")

    # population -> population (VCF )
    pop_indices = {}
    missing_samples_log = os.path.join(args.out_dir, "missing_samples_in_vcf.log")
    with open(missing_samples_log, "w") as mf:
        for pop, samples in pop2samples.items():
            idxs = []
            for s in samples:
                if s in sample_col_idx:
                    idxs.append(sample_col_idx[s])
                else:
                    mf.write(f"{pop}\t{s}\tNOT_IN_VCF\n")
            pop_indices[pop] = idxs
            print(f" {pop}: VCF {len(idxs)}/{len(samples)}", file=sys.stderr)

    types_wanted = {"DUP", "CNV"}
    if args.include_del:
        types_wanted.add("DEL")

    out_files = {}
    for pop in pop2samples:
        for t in (["DUP", "DEL"] if args.include_del else ["DUP"]):
            path = os.path.join(args.out_dir, f"{pop}.1kgp_sv.{t}.freq.bed")
            out_files[(pop, t)] = open(path, "w")

    n_records = 0
    n_used = 0
    for line in stream:
        if not line.strip() or line.startswith("#"):
            continue
        n_records += 1
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 10:
            continue
        chrom, pos, _id, _ref, _alt, _qual, filt, info_str, format_str = parts[:9]
        if args.require_pass and filt not in ("PASS", "."):
            continue

        info = parse_info(info_str)
        svtype = info.get("SVTYPE")
        if svtype not in types_wanted:
            continue

        try:
            start = int(pos)
        except ValueError:
            continue
        end = None
        if "END" in info:
            try:
                end = int(info["END"])
            except ValueError:
                end = None
        if end is None and "SVLEN" in info:
            try:
                end = start + abs(int(info["SVLEN"]))
            except ValueError:
                end = start + 1
        if end is None or end <= start:
            continue

        fmt_fields = format_str.split(":")
        gt_idx = fmt_fields.index("GT") if "GT" in fmt_fields else None
        cn_idx = (fmt_fields.index("CN") if "CN" in fmt_fields
                  else (fmt_fields.index("RD_CN") if "RD_CN" in fmt_fields else None))

        # svtype (DUP/DEL) 
        target_types = []
        if svtype == "DUP":
            target_types = [("DUP", "GT_NONREF")]
        elif svtype == "DEL" and args.include_del:
            target_types = [("DEL", "GT_NONREF")]
        elif svtype == "CNV":
            target_types = [("DUP", "CN_GT2")]
            if args.include_del:
                target_types.append(("DEL", "CN_LT2"))

        if not target_types:
            continue

        n_used += 1
        chrom_out = chrom if chrom.startswith("chr") else "chr" + chrom

        for out_type, rule in target_types:
            for pop, idxs in pop_indices.items():
                if not idxs:
                    continue
                n_carrier = 0
                n_valid = 0
                for col in idxs:
                    if col >= len(parts):
                        continue
                    sample_field = parts[col]
                    sfields = sample_field.split(":")

                    if rule == "GT_NONREF":
                        if gt_idx is None or gt_idx >= len(sfields):
                            continue
                        gt = sfields[gt_idx]
                        if gt in (".", "./.", ".|."):
                            continue
                        n_valid += 1
                        alleles = gt.replace("|", "/").split("/")
                        if any(a not in ("0", ".") for a in alleles):
                            n_carrier += 1

                    elif rule in ("CN_GT2", "CN_LT2"):
                        if cn_idx is None or cn_idx >= len(sfields):
                            continue
                        cn_val = sfields[cn_idx]
                        if cn_val in (".", ""):
                            continue
                        try:
                            cn = int(cn_val)
                        except ValueError:
                            continue
                        n_valid += 1
                        if rule == "CN_GT2" and cn > 2:
                            n_carrier += 1
                        elif rule == "CN_LT2" and cn < 2:
                            n_carrier += 1

                if n_valid == 0:
                    continue
                freq = n_carrier / n_valid
                out_files[(pop, out_type)].write(
                    f"{chrom_out}\t{start}\t{end}\t{out_type}\t{freq:.5f}\n"
                )

        if n_records % 50000 == 0:
            print(f" {n_records} , {n_used} ", file=sys.stderr)

    if proc:
        proc.wait()
    for f in out_files.values():
        f.close()

    print(f"\n. {n_records} DUP/CNV(+DEL) {n_used} ")
    print(f" : {args.out_dir}")
    print(f"VCF : {missing_samples_log}")

if __name__ == "__main__":
    main()
