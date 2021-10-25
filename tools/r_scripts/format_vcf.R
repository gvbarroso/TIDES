# Date created: 06/12/2019
# Author: Gustavo V. Barroso
# This script converts chr labels (it is part of the filter_vcfs.sh pipeline)

vcf_moms <- as.data.frame(read.table("chr_pos_moms.txt", sep = "\t"))
vcf_dads <- as.data.frame(read.table("chr_pos_dads.txt", sep = "\t"))

names(vcf_moms) <- c("CHROM", "POS")
names(vcf_dads) <- c("CHROM", "POS")

# output from write_slim_script.R
chr_limits <- read.table("chr_limits.txt", sep = "\t", header = T)

vcf_moms$CHROM <- paste("chr", findInterval(vcf_moms$POS, c(0, chr_limits$ends)), sep = "")
vcf_dads$CHROM <- paste("chr", findInterval(vcf_dads$POS, c(0, chr_limits$ends)), sep = "")

# writes files used in next step of filter_vcfs.sh pipeline
write.table(vcf_dads, "new_chr_pos_dads.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(vcf_moms, "new_chr_pos_moms.txt", col.names = F, row.names = F, quote = F, sep = "\t")
