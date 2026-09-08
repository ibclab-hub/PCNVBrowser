# Population CNV Overlap Matrix Generator
This pipeline creates continent-level summary matrices showing how many populations in each continental group have CNVs overlapping user-defined regions.

###Command Example
```
# Basic usage:
./main.sh {rawfiles_path} {query_bedfile} {outfile_path}

# With bin size:
./main.sh {rawfiles_path} {query_bedfile} {outfile_path} --bin_size 50000
```
### Input

- rawfiles_path: A directory containing CNV BED files, each with following 5 columns(chr, start, end, type, frequency)
    You may also use **Combined CNV BED files** exported from [PCNVBrowser](http://biweb.konkuk.ac.kr/PCNVBrowser/Data_download/index.html).
- query_bedfile: A BED file containing CNV regions of interest with 4 columns(chr, start, end, type)
- outfile_path: Output directory for continent-wise matrices
- --bin_size (optional): Bin size (bp) used to subdivide each CNV region (default: 50000 bp)

### Output

Matrix files are generated only for chromosomes present in the query BED file (one file per chromosome).
Example output filenames:

```
  continent_matrix.chr1.tsv
  continent_matrix.chr2.tsv
  ...
```

Each row represents a 50 kb bin:
| Region             | AFR | EUR | EAS | SAS | AMR |
|--------------------|-----|-----|-----|-----|-----|
| chr1:100000-150000 | 2/7 | 1/5 | 0/5 | 0/5 | 1/4 |


- Each cell shows the number of populations (out of N) within each continent that have a CNV overlapping the bin.
- `/N` denotes the total number of populations considered per continent.
- Bins with no overlaps in any continent are omitted. If none of the query regions overlap population CNVs, the output may be empty.

### Dependencies

- Python 3
- bedtools
