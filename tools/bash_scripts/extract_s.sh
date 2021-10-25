#!/bin/bash

# this script starts from the output files of filter_vcfs.sh

# filters out 2nd instance of multiallelic sites (since they have different s values)
zcat fathers.vcf.gz | grep -v 'MULTIALLELIC' > fathers_single_allelic.vcf 
zcat mothers.vcf.gz | grep -v 'MULTIALLELIC' > mothers_single_allelic.vcf

# selects non-synonymous mutations
grep 'MT=1' fathers_single_allelic.vcf > fathers_single_allelic_mt1.vcf 
grep 'MT=1' mothers_single_allelic.vcf > mothers_single_allelic_mt1.vcf

cat moms_header.txt mothers_single_allelic_mt1.vcf | bgzip > mothers_single_allelic_mt1.vcf.gz 
cat dads_header.txt fathers_single_allelic_mt1.vcf | bgzip > fathers_single_allelic_mt1.vcf.gz

# selects synonymous mutations
grep 'MT=2' fathers_single_allelic.vcf > fathers_single_allelic_mt2.vcf 
grep 'MT=2' mothers_single_allelic.vcf > mothers_single_allelic_mt2.vcf

cat moms_header.txt mothers_single_allelic_mt2.vcf | bgzip > mothers_single_allelic_mt2.vcf.gz 
cat dads_header.txt fathers_single_allelic_mt2.vcf | bgzip > fathers_single_allelic_mt2.vcf.gz

# selects CHROM and POS columns for NS mutations, repeating POS column to create BEDGRAPH file
bcftools view -h mothers_single_allelic_mt1.vcf.gz > moms_header_mt1.txt 
bcftools view -h mothers_single_allelic_mt1.vcf.gz > dads_header_mt1.txt

header_size_moms=$(wc -l moms_header_mt1.txt | cut -d ' ' -f 1)
cut_size_moms=$(($header_size_moms))

header_size_dads=$(wc -l dads_header_mt1.txt | cut -d ' ' -f 1)
cut_size_dads=$(($header_size_dads)) 

zcat fathers_single_allelic_mt1.vcf.gz | tail -n +"$cut_size_dads" | cut -f 1-2 > f12.txt 
zcat mothers_single_allelic_mt1.vcf.gz | tail -n +"$cut_size_moms" | cut -f 1-2 > m12.txt

cut -f 2 f12.txt > f2.txt 
cut -f 2 m12.txt > m2.txt
paste f12.txt f2.txt  > f_coords.txt 
paste m12.txt m2.txt > m_coords.txt

# selects s values for mutations
zcat fathers_single_allelic_mt1.vcf.gz | tail -n +"$cut_size_dads" |  cut -f8 | grep -oE $'[;\t]S=[^;\t]+' | awk 'BEGIN { FS = "="; } { print $2 }' > s_fathers.txt 
zcat mothers_single_allelic_mt1.vcf.gz | tail -n +"$cut_size_moms" | cut -f8 | grep -oE $'[;\t]S=[^;\t]+' | awk 'BEGIN { FS = "="; } { print $2 }' > s_mothers.txt

# assembles BEDGRAPH files with s values 
paste f_coords.txt s_fathers.txt > dfe_fathers.bedgraph 
paste m_coords.txt s_mothers.txt  > dfe_mothers.bedgraph

# to help compute the (weighted) mean s of segregating mutations
zcat fathers_single_allelic_mt1.vcf.gz | tail -n +"$cut_size_dads" | cut -f8 | grep -oE $'[;\t]AC=[^;\t]+' | awk 'BEGIN { FS = "="; } { print $2 }' > ac_fathers.txt 
zcat mothers_single_allelic_mt1.vcf.gz | tail -n +"$cut_size_moms" | cut -f8 | grep -oE $'[;\t]AC=[^;\t]+' | awk 'BEGIN { FS = "="; } { print $2 }' > ac_mothers.txt

paste dfe_fathers.bedgraph ac_fathers.txt | cut -f 1,2,4,5 > dfe_ac_fathers.tsv
paste dfe_mothers.bedgraph ac_mothers.txt | cut -f 1,2,4,5 > dfe_ac_mothers.tsv

rm *.vcf
