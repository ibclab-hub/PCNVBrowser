# Role: Aggregates the per-population overlap matrix into five continental super-population summaries (Application 3, Table 2).
import pandas as pd

AFR = ["ACB","ASW","ESN","GWD","LWK","MSL","YRI"] # 7
AMR = ["CLM","MXL","PEL","PUR"] # 4
EAS = ["CDX","CHB","CHS","JPT","KHV"] # 5
EUR = ["CEU","FIN","GBR","IBS","TSI"] # 5
SAS = ["BEB","GIH","ITU","PJL","STU"] # 5

groups = {"East Asian": EAS, "South Asian": SAS, "European": EUR,
          "African": AFR, "American": AMR}

for tool_file in ["CN.MOPS_25K_population_matrix.tsv",
                  "CNVKIT_5K_population_matrix.tsv",
                  "FREEC_50K_population_matrix.tsv",
                  "READDEPTH_population_matrix.tsv"]:
    df = pd.read_csv(tool_file, sep="\t", index_col="Region")
    n = len(df)
    print(tool_file)
    for cont, pops in groups.items():
        pops_present = [p for p in pops if p in df.columns]
        overlap_count = (df[pops_present].sum(axis=1) > 0).sum()
        print(f"  {cont}: {overlap_count}/{n} ({overlap_count/n*100:.2f}%)")
