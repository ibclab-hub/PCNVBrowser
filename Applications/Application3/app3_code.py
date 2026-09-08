# Role: Builds a presence/absence matrix between the 15 Korean-specific CNVs and each of the 26 population CNV tracks (Application 3).
import pandas as pd
import sys
import glob
import os

def load_query_cnvs(bed_file):
    queries = []
    with open(bed_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            cnv_type = parts[3] if len(parts) > 3 else "."
            queries.append((chrom, start, end, cnv_type))

    print(f"[DEBUG] Loaded {len(queries)} query CNVs from {bed_file}")
    for q in queries[:5]:
        print(f"  {q[0]}:{q[1]}-{q[2]} ({q[3]})")
    if len(queries) > 5:
        print(f"  ... and {len(queries)-5} more")

    return queries

def load_population_beds(bed_dir):
    pattern = os.path.join(bed_dir, "*.all_0.05.bed")
    files = sorted(glob.glob(pattern))

    pop_dict = {}
    for f in files:
        pop_name = os.path.basename(f).split(".")[0]  # ACB.all_0.05.bed -> ACB
        pop_dict[pop_name] = f

    print(f"[DEBUG] Found {len(pop_dict)} population BED files in {bed_dir}")
    for pop, path in pop_dict.items():
        print(f"  {pop}: {path}")

    return pop_dict

def load_bed_regions(bed_file):
    regions = []
    with open(bed_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            regions.append((parts[0], int(parts[1]), int(parts[2])))
    return regions

def has_overlap(start_a, end_a, start_b, end_b):
    return min(end_a, end_b) > max(start_a, start_b)

def make_population_matrix(pop_dict, queries, output_file="population_matrix.tsv"):
    row_labels = [f"{chrom}:{start}-{end}_{cnv_type}" for chrom, start, end, cnv_type in queries]
    matrix = pd.DataFrame(0, index=row_labels, columns=pop_dict.keys())

    for pop, bed_file in pop_dict.items():
        subjects = load_bed_regions(bed_file)
        print(f"[DEBUG] {pop}: loaded {len(subjects)} entries")

        match_count = 0
        for i, (chrom_q, start_q, end_q, _) in enumerate(queries):
            matched = False
            for chrom_b, start_b, end_b in subjects:
                if chrom_b != chrom_q:
                    continue
                if has_overlap(start_q, end_q, start_b, end_b):
                    matched = True
                    break # 1/0 
            if matched:
                matrix.at[row_labels[i], pop] = 1
                match_count += 1
        print(f"[DEBUG] {pop}: {match_count}/{len(queries)} queries matched")

    matrix.index.name = "Region"
    matrix.to_csv(output_file, sep="\t")
    print(f"[✓] Matrix saved to: {output_file}")
    return matrix

# main
query_file = sys.argv[1]
bed_dir = sys.argv[2]
output_file = sys.argv[3]

queries = load_query_cnvs(query_file)
pop_dict = load_population_beds(bed_dir)
make_population_matrix(pop_dict, queries, output_file=output_file)