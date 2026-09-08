#!/usr/bin/python3
import os
import argparse
import pandas as pd
from collections import defaultdict

def create_50k_bins(start, end, bin_size=50000):
    bin_start = (start // bin_size) * bin_size
    bin_end = ((end // bin_size) + 1) * bin_size
    return [(b, b + bin_size) for b in range(bin_start, bin_end, bin_size)]

def make_suffix(matrix):
    suffix_map = {
        'EAS': '/5',
        'EUR': '/5',
        'SAS': '/5',
        'AMR': '/4',
        'AFR': '/7'
    }
    for continent, suffix in suffix_map.items():
        matrix[continent] = matrix[continent].apply(lambda x: f"{x}{suffix}" if x != 0 else "0")
    matrix.index.name = "Region"
    
    return matrix

def make_continent_matrix(bed_dict, region, bin_size):
    chrom, coord = region.split(":")
    start_region, end_region = map(int, coord.split("-"))

    bins = create_50k_bins(start_region, end_region, bin_size)
    bin_labels = [f"{chrom}:{start}-{end}" for start, end in bins]
    matrix = pd.DataFrame(0, index=bin_labels, columns=bed_dict.keys())

    for continent, bed_file in bed_dict.items():
        with open(bed_file) as f:
            for line in f:
                chrom_b, start_b, end_b, count = line.strip().split("\t")
                if chrom_b != chrom:
                    continue
                start_b = int(start_b)
                end_b = int(end_b)
                count = int(count)

                for i, (bin_start, bin_end) in enumerate(bins):
                    if start_b < bin_end and end_b > bin_start:
                        matrix.iloc[i][continent] = max(matrix.iloc[i][continent], count)

    matrix = matrix[(matrix.T != 0).any()]
 
    return matrix

def get_chrom_regions(bed_file):
    chrom_ranges = defaultdict(list)

    with open(bed_file) as f:
        for line in f:
            chrom, start, end, _ = line.strip().split("\t")
            start = int(start)
            end = int(end)
            
            chrom_ranges[chrom].append((start, end))

    return chrom_ranges

def main():
    parser = argparse.ArgumentParser(description="Generate CNV overlap matrix by continent for given CNV regions.")
    parser.add_argument("bed_file", help="BED file with regions of interest (e.g., Korean CNVs)")
    parser.add_argument("multi_dir", help="Directory containing multiinter continent BEDs")
    parser.add_argument("out_dir", help="Directory to save output matrices")
    parser.add_argument("--bin_size", type=int, default=50000, help="Bin size for region subdivision (default: 50000)")

    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    bed_dict = {
        'AFR': os.path.join(args.multi_dir, 'AFR.multi.bed'),
        'AMR': os.path.join(args.multi_dir, 'AMR.multi.bed'),
        'EAS': os.path.join(args.multi_dir, 'EAS.multi.bed'),
        'EUR': os.path.join(args.multi_dir, 'EUR.multi.bed'),
        'SAS': os.path.join(args.multi_dir, 'SAS.multi.bed'),
    }

    chrom_ranges = get_chrom_regions(args.bed_file)

    for chrom, region_list in chrom_ranges.items():
        result_matrix=pd.DataFrame()
        for start, end in region_list:
            region = f"{chrom}:{start}-{end}"
            tmp = make_continent_matrix(bed_dict, region, args.bin_size)
            
            if result_matrix.empty:
                result_matrix = tmp
            else:
                result_matrix = pd.concat([result_matrix, tmp])
                result_matrix = result_matrix.groupby(result_matrix.index).max()
        
        result_matrix = make_suffix(result_matrix)
        output_file = os.path.join(args.out_dir, f"continent_matrix.{chrom}.tsv")
        
        result_matrix.to_csv(output_file, sep="\t")
        print(f"-Matrix saved to: {output_file}")

if __name__ == "__main__":
    main()

