#!/usr/bin/env python3
"""
recall_by_truth_size.py

DEL recall이 낮게 나온 게 "caller 성능이 나빠서"인지 "truth set에 caller가
애초에 탐지 못하는 초소형 변이가 압도적으로 많아서(size mismatch)"인지
구분하기 위해, truth DEL을 크기 구간별로 쪼개서 recall을 따로 계산.

구간(기본): 50-500bp, 500bp-1kb, 1-5kb, 5-25kb, >=25kb
(PCNVBrowser 최소 bin이 5000bp인 걸 감안해 5kb를 기준선 중 하나로 포함)

사용 예:
  python3 recall_by_truth_size.py \
      --truth-dir truth_beds \
      --data-dir /var/www/html/PCNVBrowser/data/samples \
      --samples HG00438,HG00621,...,NA12878 \
      --resolution 25000 \
      --tool READDEPTH \
      --out-dir ./size_stratified_results
"""
import argparse
import csv
import os
import subprocess
import sys
import tempfile

SIZE_BINS = [(50, 500), (500, 1000), (1000, 5000), (5000, 25000), (25000, float("inf"))]


def bin_label(lo, hi):
    if hi == float("inf"):
        return f">={lo}bp"
    return f"{lo}-{hi}bp"


def find_raw_bed(data_dir, sample, tool, resolution):
    sample_dir = os.path.join(data_dir, sample)
    if tool == "READDEPTH":
        path = os.path.join(sample_dir, f"READDEPTH.{sample}.sorted.bed")
    else:
        path = os.path.join(sample_dir, f"{tool}.{resolution}.{sample}.sorted.bed")
    return path if os.path.isfile(path) else None


def extract_del_sorted(src_path, dst_path):
    rows = []
    with open(src_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 4 or parts[3] != "DEL":
                continue
            rows.append((parts[0], int(parts[1]), int(parts[2])))
    rows.sort()
    with open(dst_path, "w") as out:
        for chrom, start, end in rows:
            out.write(f"{chrom}\t{start}\t{end}\n")
    return len(rows)


def split_truth_by_size(truth_del_path, dst_dir, sample):
    """truth DEL.bed를 크기 구간별로 쪼개서 저장, {bin_label: (path, n)} 반환."""
    bins = {b: [] for b in SIZE_BINS}
    with open(truth_del_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            start, end = int(parts[1]), int(parts[2])
            size = end - start
            for lo, hi in SIZE_BINS:
                if lo <= size < hi:
                    bins[(lo, hi)].append((parts[0], start, end))
                    break
    out = {}
    for (lo, hi), rows in bins.items():
        rows.sort()
        path = os.path.join(dst_dir, f"{sample}_truth_{lo}_{hi}.bed")
        with open(path, "w") as f:
            for chrom, s, e in rows:
                f.write(f"{chrom}\t{s}\t{e}\n")
        out[(lo, hi)] = (path, len(rows))
    return out


def split_calls_by_size(call_bed_path, dst_dir, sample, tag):
    """우리 콜(call_bed_path)을 자기 자신의 크기로 구간별로 쪼갬 -> precision 계산용."""
    bins = {b: [] for b in SIZE_BINS}
    with open(call_bed_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            start, end = int(parts[1]), int(parts[2])
            size = end - start
            for lo, hi in SIZE_BINS:
                if lo <= size < hi:
                    bins[(lo, hi)].append((parts[0], start, end))
                    break
    out = {}
    for (lo, hi), rows in bins.items():
        rows.sort()
        path = os.path.join(dst_dir, f"{sample}_{tag}_calls_{lo}_{hi}.bed")
        with open(path, "w") as f:
            for chrom, s, e in rows:
                f.write(f"{chrom}\t{s}\t{e}\n")
        out[(lo, hi)] = (path, len(rows))
    return out


def reciprocal_frac(path_a, path_b, n_a, frac=0.5):
    if n_a == 0:
        return None
    try:
        result = subprocess.run(
            ["bedtools", "intersect", "-u", "-f", str(frac), "-r", "-a", path_a, "-b", path_b],
            capture_output=True, text=True, check=True,
        )
    except subprocess.CalledProcessError:
        return None
    n_hit = len([l for l in result.stdout.splitlines() if l.strip()])
    return n_hit / n_a


def restrict_to_confident(call_bed, confident_bed, dst_path):
    """call_bed 중 confident_bed와 조금이라도 겹치는 콜만 남김.
    Confident 밖 영역은 애초에 truth가 없어 판단 불가능하므로, recall/precision
    계산에서 부당하게 false positive로 잡히지 않도록 사전 제외.
    (validate_against_truth.py와 동일한 로직 재사용 — S5/S6 일관성 유지)"""
    if confident_bed is None or not os.path.isfile(confident_bed):
        return call_bed
    try:
        result = subprocess.run(
            ["bedtools", "intersect", "-u", "-a", call_bed, "-b", confident_bed],
            capture_output=True, text=True, check=True,
        )
    except subprocess.CalledProcessError:
        return call_bed
    lines = [l for l in result.stdout.splitlines() if l.strip()]
    if not lines:
        return None
    with open(dst_path, "w") as f:
        for l in lines:
            f.write(l + "\n")
    return dst_path


def find_confident_bed(truth_dir, sample):
    path = os.path.join(truth_dir, sample, "CONFIDENT.bed")
    return path if os.path.isfile(path) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth-dir", required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--resolution", default="25000")
    ap.add_argument("--tool", default="READDEPTH",
                     help="어느 tool의 raw call로 recall 계산할지 (기본 READDEPTH, "
                          "다른 tool 이름 또는 'ALL_4TOOL_UNION' 지정 가능)")
    ap.add_argument("--reciprocal-frac", type=float, default=0.5)
    ap.add_argument("--out-dir", default="./size_stratified_results")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    samples = [s.strip() for s in args.samples.split(",") if s.strip()]

    out_path = os.path.join(args.out_dir, f"recall_by_size_{args.tool}.csv")
    precision_out_path = os.path.join(args.out_dir, f"precision_by_size_{args.tool}.csv")

    with open(out_path, "w", newline="") as f, \
         open(precision_out_path, "w", newline="") as fp, \
         tempfile.TemporaryDirectory() as tmpdir:

        writer = csv.writer(f)
        writer.writerow(["sample", "size_bin", "n_truth_in_bin", "n_calls", "recall"])
        writer_p = csv.writer(fp)
        writer_p.writerow(["sample", "size_bin", "n_calls_in_bin", "n_truth_total", "precision"])

        for sample in samples:
            truth_del = os.path.join(args.truth_dir, sample, "DEL.bed")
            if not os.path.isfile(truth_del):
                continue

            size_bins = split_truth_by_size(truth_del, tmpdir, sample)
            n_truth_total = sum(n for _, n in size_bins.values())

            if args.tool == "ALL_4TOOL_UNION":
                tools = ["CN.MOPS", "CNVKIT", "FREEC", "READDEPTH"]
                paths = []
                for t in tools:
                    raw = find_raw_bed(args.data_dir, sample, t, args.resolution)
                    if raw:
                        p = os.path.join(tmpdir, f"{sample}_{t}_DEL.bed")
                        extract_del_sorted(raw, p)
                        paths.append(p)
                combined = os.path.join(tmpdir, f"{sample}_ALL_DEL.bed")
                rows = set()
                for p in paths:
                    with open(p) as pf:
                        for l in pf:
                            if l.strip():
                                rows.add(l.strip())
                with open(combined, "w") as cf:
                    for r in sorted(rows):
                        cf.write(r + "\n")
                call_bed = combined
                n_calls = len(rows)
            else:
                raw = find_raw_bed(args.data_dir, sample, args.tool, args.resolution)
                if raw is None:
                    continue
                call_bed = os.path.join(tmpdir, f"{sample}_{args.tool}_DEL.bed")
                n_calls = extract_del_sorted(raw, call_bed)

            # confident region 제한 (validate_against_truth.py와 동일하게 적용):
            # confident 밖 콜은 애초에 판단 불가능하므로 recall/precision 계산 전에 제외.
            confident_bed = find_confident_bed(args.truth_dir, sample)
            if confident_bed and n_calls > 0:
                restricted = restrict_to_confident(
                    call_bed, confident_bed,
                    os.path.join(tmpdir, f"{sample}_{args.tool}_DEL_conf.bed")
                )
                if restricted is None:
                    n_calls = 0
                else:
                    call_bed = restricted
                    with open(call_bed) as _f:
                        n_calls = sum(1 for _ in _f)

            # (A) recall: truth를 크기별로 쪼개서, 우리 콜 전체(call_bed) 대비
            for (lo, hi), (truth_bin_path, n_truth_bin) in size_bins.items():
                label = bin_label(lo, hi)
                if n_truth_bin == 0 or n_calls == 0:
                    writer.writerow([sample, label, n_truth_bin, n_calls, "NA"])
                    continue
                recall = reciprocal_frac(truth_bin_path, call_bed, n_truth_bin, args.reciprocal_frac)
                writer.writerow([sample, label, n_truth_bin,
                                  n_calls, f"{recall:.4f}" if recall is not None else "NA"])

            # (B) precision: 우리 콜을 크기별로 쪼개서, truth 전체(truth_del) 대비
            if n_calls > 0:
                call_bins = split_calls_by_size(call_bed, tmpdir, sample, args.tool)
                for (lo, hi), (call_bin_path, n_calls_bin) in call_bins.items():
                    label = bin_label(lo, hi)
                    if n_calls_bin == 0 or n_truth_total == 0:
                        writer_p.writerow([sample, label, n_calls_bin, n_truth_total, "NA"])
                        continue
                    precision = reciprocal_frac(call_bin_path, truth_del, n_calls_bin,
                                                 args.reciprocal_frac)
                    writer_p.writerow([sample, label, n_calls_bin, n_truth_total,
                                        f"{precision:.4f}" if precision is not None else "NA"])

            print(f"{sample} 완료", file=sys.stderr)

    print(f"\nrecall 결과: {out_path}")
    print(f"precision 결과: {precision_out_path}")


if __name__ == "__main__":
    main()
