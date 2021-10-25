#!/bin/bash

# saves the headers
bcftools view -h potential_mothers.vcf.gz > moms_header.txt 
bcftools view -h potential_fathers.vcf.gz > dads_header.txt

header_size_moms=$(wc -l moms_header.txt | cut -d ' ' -f 1)
cut_size_moms=$(($header_size_moms - 2)) # -2 since bcftools adds 3 lines to header

header_size_dads=$(wc -l dads_header.txt | cut -d ' ' -f 1)
cut_size_dads=$(($header_size_dads - 2)) # -2 since bcftools adds 3 lines to header

# cuts headers off
zcat potential_mothers.vcf.gz | tail -n +"$cut_size_moms" > potential_mothers_formatted.vcf 
zcat potential_fathers.vcf.gz | tail -n +"$cut_size_dads" > potential_fathers_formatted.vcf

# grabs CHROM and POS columns
cut -f 1-2 potential_mothers_formatted.vcf > chr_pos_moms.txt 
cut -f 1-2 potential_fathers_formatted.vcf > chr_pos_dads.txt

# formats CHROM column
# (this script reads 'chr_limits.txt', 'chr_pos_moms.txt' and 'chr_pos_dads.txt')
# (it outputs 'new_chr_pos_moms.txt' and 'new_chr_pos_dads.txt'
Rscript format_vcf.R

# edits CHROM column inside 'potential_mothers_formatted.vcf' and 'potential_fathers_formatted.vcf'
awk 'FNR==NR{a[NR]=$1;next}{OFS="\t";$1=a[FNR]}1' new_chr_pos_moms.txt potential_mothers_formatted.vcf > moms.vcf 
awk 'FNR==NR{a[NR]=$1;next}{OFS="\t";$1=a[FNR]}1' new_chr_pos_dads.txt potential_fathers_formatted.vcf > dads.vcf

# assembles formatted VCFs with original headers
cat moms_header.txt moms.vcf | bgzip > mothers.vcf.gz 
cat dads_header.txt dads.vcf | bgzip > fathers.vcf.gz

rm *.vcf
