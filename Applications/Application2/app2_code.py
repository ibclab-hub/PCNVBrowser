#!/usr/bin/python
# Role: Filters gnomAD-SV common CNVs and reported PD-associated CNVs against European population CNV tracks (Application 2, Table 1).
import sys
import pandas as pd

def compare_all(a_df, b_df, out="out.csv"):
    with open(out, "w") as fw:
        for _, b_row in b_df.iterrows():
            b_chrom = b_row['chrom']
            b_start = int(b_row['start'])
            b_end = int(b_row['end'])

            overlapping = a_df[
                (a_df["chrom"] == b_chrom) &
                (a_df["start"] < b_end) &
                (a_df["end"] > b_start)
            ]

            if len(overlapping) == 0:
                fw.write(f"{b_chrom}\t{b_start}\t{b_end}\n")

def run(b_file):
    b_df = pd.read_csv(b_file, sep='\t', header=None,
                        names=['chrom', 'start', 'end', 'type'])
    eur_dir = "./data/EUR_data"

    for tool in ['READDEPTH', 'CN.MOPS_25000', 'CNVKIT_5000', 'FREEC_50000']:
        a_df = pd.read_csv(f"{eur_dir}/EUR.{tool}.ALL.0.05.bed",
                           sep='\t', header=None,
                           names=['chrom', 'start', 'end'])

        out_file = f"./tmp/intersect.EUR.{tool}.ALL.csv"
        compare_all(a_df, b_df, out_file)

        result = pd.read_csv(out_file, sep='\t', header=None,
                             names=['chrom', 'start', 'end'])
        print(tool, len(result), len(result) / len(b_df) * 100)

if __name__ == '__main__':
    run(sys.argv[1])