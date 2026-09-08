#!/usr/bin/env python3
# Role: Builds standardized truth BED files from HPRC dipcall assemblies (26 additional samples).

import argparse
import bisect
import gzip
import os
from collections import defaultdict

def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)

def normalize_chrom(chrom):
    if chrom.startswith("chr"):
        return chrom
    return "chr" + chrom

def load_confident_regions(path):
    regions = defaultdict(list)

    with open_maybe_gzip(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            p = line.rstrip().split()
            if len(p) < 3:
                continue

            chrom = normalize_chrom(p[0])
            start = int(p[1])
            end = int(p[2])

            regions[chrom].append((start, end))

    for chrom in regions:
        regions[chrom].sort()

    return regions

def build_region_index(regions):
    starts = {}

    for chrom, arr in regions.items():
        starts[chrom] = [x[0] for x in arr]

    return starts

def fully_in_confident_region(chrom, start, end, regions, starts):
    if chrom not in regions:
        return False

    arr = regions[chrom]
    pos = starts[chrom]

    i = bisect.bisect_right(pos, start) - 1

    if i < 0:
        return False

    r_start, r_end = arr[i]

    return r_start <= start and end <= r_end

def write_bed(rows, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)

    rows = sorted(set(rows))

    with open(path, "w") as f:
        for chrom, start, end in rows:
            f.write(f"{chrom}\t{start}\t{end}\n")

    return len(rows)

def main():

    ap = argparse.ArgumentParser()

    ap.add_argument("--sample", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--confident-bed", required=True)

    ap.add_argument(
        "--min-size",
        type=int,
        default=50,
        help="minimum absolute indel size (default: 50 bp)"
    )

    ap.add_argument(
        "--out-dir",
        default="truth_beds"
    )

    ap.add_argument(
        "--include-filtered",
        action="store_true",
        help=(
            "FILTER != PASS/. variant . "
            " PASS '.' ."
        )
    )

    args = ap.parse_args()

    # ----------------------------------------------------------
    # confident regions
    # ----------------------------------------------------------

    regions = load_confident_regions(args.confident_bed)
    starts = build_region_index(regions)

    rows_del = []
    rows_ins = []

    stats = defaultdict(int)

    # ----------------------------------------------------------
    # VCF
    # ----------------------------------------------------------

    with open_maybe_gzip(args.vcf) as f:

        for line in f:

            if not line.strip() or line.startswith("#"):
                continue

            stats["vcf_records"] += 1

            p = line.rstrip().split("\t")

            if len(p) < 8:
                stats["malformed"] += 1
                continue

            chrom = normalize_chrom(p[0])

            try:
                pos = int(p[1])
            except ValueError:
                stats["bad_pos"] += 1
                continue

            ref = p[3]
            alt_field = p[4]
            filt = p[6]

            # --------------------------------------------------
            # dipcall filtered records
            # --------------------------------------------------

            if not args.include_filtered:
                if filt not in ("PASS", "."):
                    stats["filtered_vcf"] += 1
                    continue

            # --------------------------------------------------
            # multi-allelic ALT 
            # --------------------------------------------------

            for alt in alt_field.split(","):

                # symbolic/breakend sequence parser 
                if (
                    alt == "."
                    or alt.startswith("<")
                    or "[" in alt
                    or "]" in alt
                    or "*" == alt
                ):
                    stats["symbolic_or_bnd"] += 1
                    continue

                diff = len(alt) - len(ref)

                # ==================================================
                # DELETION
                #
                # VCF:
                # POS=100 REF=ATTTT ALT=A
                #
                # POS is 1-based anchor.
                # sequence = reference 101-104
                #
                # BED:
                # start=100, end=104
                # ==================================================

                if diff <= -args.min_size:

                    del_len = -diff

                    bed_start = pos
                    bed_end = pos + del_len

                    stats["del_size_pass"] += 1

                    if fully_in_confident_region(
                        chrom,
                        bed_start,
                        bed_end,
                        regions,
                        starts,
                    ):
                        rows_del.append(
                            (chrom, bed_start, bed_end)
                        )
                        stats["del_confident"] += 1
                    else:
                        stats["del_outside_confident"] += 1

                # ==================================================
                # INSERTION
                #
                # insertion reference interval .
                # marker 1 bp BED .
                #
                # DUP truth .
                # ==================================================

                elif diff >= args.min_size:

                    ins_len = diff

                    # anchor base insertion
                    bed_start = pos
                    bed_end = pos + 1

                    stats["ins_size_pass"] += 1

                    if fully_in_confident_region(
                        chrom,
                        bed_start,
                        bed_end,
                        regions,
                        starts,
                    ):
                        rows_ins.append(
                            (chrom, bed_start, bed_end)
                        )
                        stats["ins_confident"] += 1
                    else:
                        stats["ins_outside_confident"] += 1

                else:
                    stats["small_variant"] += 1

    # ----------------------------------------------------------
    # write
    # ----------------------------------------------------------

    sample_dir = os.path.join(args.out_dir, args.sample)

    del_path = os.path.join(sample_dir, "DEL.bed")
    ins_path = os.path.join(sample_dir, "INS.bed")

    nd = write_bed(rows_del, del_path)
    ni = write_bed(rows_ins, ins_path)

    # confident region / 
    confident_out = os.path.join(
        sample_dir,
        "CONFIDENT.bed"
    )

    confident_rows = []

    for chrom, arr in regions.items():
        for start, end in arr:
            confident_rows.append(
                (chrom, start, end)
            )

    nc = write_bed(
        confident_rows,
        confident_out
    )

    print()
    print(f"===== {args.sample} =====")
    print(f"VCF records              : {stats['vcf_records']}")
    print(f"VCF filtered             : {stats['filtered_vcf']}")
    print(f"small variants           : {stats['small_variant']}")
    print()
    print(f"DEL >= {args.min_size} bp")
    print(f"  size pass              : {stats['del_size_pass']}")
    print(f"  confident              : {stats['del_confident']}")
    print(f"  outside confident      : {stats['del_outside_confident']}")
    print()
    print(f"INS >= {args.min_size} bp")
    print(f"  size pass              : {stats['ins_size_pass']}")
    print(f"  confident              : {stats['ins_confident']}")
    print(f"  outside confident      : {stats['ins_outside_confident']}")
    print()
    print(f"DEL truth                : {nd} -> {del_path}")
    print(f"INS reference-only       : {ni} -> {ins_path}")
    print(f"confident regions        : {nc} -> {confident_out}")
    print()
    print("NOTE: INS.bed must NOT be used as DUP truth.")
    print()

if __name__ == "__main__":
    main()
