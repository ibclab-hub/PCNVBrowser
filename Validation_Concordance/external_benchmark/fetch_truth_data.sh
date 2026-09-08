#!/bin/bash
# Role: Downloads the external reference/truth data sets used for validation.
# fetch_truth_data.sh
# 42 validation truth 3 
#
# 1) NA12878 (Parikh et al. 2016) - URL, DEL, GRCh37
# 2) HGSVC2 freeze3 (Ebert et al. 2021) - URL, HG00514/HG00733/NA19240 ,
# DEL+INS( DUP INS ), merged TSV 
# 3) HPRC Year1 dipcall - S3 , 

set -e
mkdir -p truth_raw
cd truth_raw

echo "===== 1) NA12878 (Parikh et al. 2016) ====="
wget -c "ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/svclassify_Manuscript/Supplementary_Information/Personalis_1000_Genomes_deduplicated_deletions.bed" \
    -O NA12878.parikh2016.deletions.GRCh37.bed
echo "-> DEL , build=GRCh37. PCNVBrowser GRCh38 liftover (lift_na12878.sh )"

echo ""
echo "===== 2) HGSVC2 freeze3 (Ebert et al. 2021) ====="
wget -c "https://zenodo.org/records/4268828/files/variants_freeze3_sv_insdel.tsv.gz?download=1" \
    -O hgsvc2_freeze3_sv_insdel.tsv.gz
echo "-> 32 TSV . :"
echo "     zcat hgsvc2_freeze3_sv_insdel.tsv.gz | head -3"
echo " ( sample genotype HG00514/HG00733/NA19240 ) "
echo " parse_hgsvc2.py SAMPLE_COL_PATTERN ."

echo ""
echo "===== 3) HPRC Year1 dipcall ====="
echo "-> S3 . :"
echo ""
echo "   aws s3 ls --no-sign-request s3://human-pangenomics/working/HPRC/HG00438/ --recursive | grep -i dip"
echo ""
echo " '*.dip.vcf.gz' '*.dip.bed' ,"
echo " download_hprc.sh DIP_VCF_PATTERN / DIP_BED_PATTERN ."
echo " (aws-cli : pip install awscli --break-system-packages , "
echo " https://s3.amazonaws.com/human-pangenomics/working/HPRC/HG00438/ )"
