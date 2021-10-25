# Date created: 02/12/2019
# Author: Gustavo V. Barroso
# Retrieves and organizes data from Human Genome, to set up SLiM simulations

library("biomaRt")
library("GenomicRanges")
library("tidyverse")
library("plyr")

setwd("~/Data/TIDES/")

#######################
#
# Loads & Filters data
#
#######################

# gets exome data from human genome (accessed on January 5 2020 to perform simulations)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # hg38

human_exome_ensembl <- getBM(mart = ensembl, filters = c("chromosome_name", "with_ccds"), values = list(c(1:22), T),
                             attributes = c("ensembl_gene_id", "ensembl_exon_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end"))

human_exome <- dplyr::arrange(.data = human_exome_ensembl, chromosome_name, exon_chrom_start)

# filters out overlapping exons / genes
human_exome <- ddply(.data = human_exome, .variables = "chromosome_name", .fun = distinct, exon_chrom_start, .keep_all = T)
human_exome <- ddply(.data = human_exome, .variables = "chromosome_name", .fun = distinct, exon_chrom_end, .keep_all = T)

# we use genomic ranges objects for speed
human_exome_gr <- makeGRangesFromDataFrame(human_exome, keep.extra.columns = T)

# iteratively gets rid of the final overlapping exons
curr_num_ranges <- length(ranges(human_exome_gr))
prev_num_ranges <- curr_num_ranges + 1
while(prev_num_ranges > curr_num_ranges) {
  human_exome_gr <- human_exome_gr[unique(findOverlaps(human_exome_gr, type = "any", select = "first"))];
  prev_num_ranges <- curr_num_ranges
  curr_num_ranges <- length(ranges(human_exome_gr))
}

# merges contiguous exons from different genes 
cont_gr <- as.data.frame(human_exome_gr)
good_rows <- rep(T, nrow(cont_gr))
for(i in 2:nrow(cont_gr)) {
  if(cont_gr[i, 2] == (cont_gr[i - 1, 3] + 1)) {
    cont_gr[i - 1, 3] <- cont_gr[i, 3]
    good_rows[i] <- F
  }
}
cont_gr <- cont_gr[good_rows, c(1:3, 6:7)]
row.names(cont_gr) <- 1:nrow(cont_gr)

# safety first
cont_gr$seqnames <- as.numeric(cont_gr$seqnames)
cont_gr$start <- as.numeric(cont_gr$start)
cont_gr$end <- as.numeric(cont_gr$end)

# to match chr labels in rec. maps
for(i in 1:22) { cont_gr$seqnames[which(cont_gr$seqnames == i)] <- paste("chr", i, sep = "") }

# updates GR
human_exome_gr <- makeGRangesFromDataFrame(cont_gr, keep.extra.columns = T)

# female and male recombination maps binned into windows ~100kb (DECODE 2019)
# raw maps downloaded from https://science.sciencemag.org/content/363/6425/eaau1043 on November 19, 2019

## female
frmap <- read.table("human_genome/female_rmap_DECODE.bedgraph", header = F, sep = "\t")
names(frmap) <-  c("chrom", "chromStart", "chromEnd", "rate")
frmap <- frmap[-which(frmap$chrom %in% "chrX"),]
frmap$chr.id <- 0
for(i in 1:22) { frmap$chr.id[which(frmap$chrom == paste("chr", i, sep = ""))] <- i }
frmap <- dplyr::arrange(.data = frmap, chr.id, chromStart)
frmap <- frmap[-which(names(frmap) == "chr.id")]
frmap$rate <- frmap$rate * 1e-8 # converts from cM/Mb to r

## male
mrmap <- read.table("human_genome/male_rmap_DECODE.bedgraph", header = F, sep = "\t")
names(mrmap) <-  c("chrom", "chromStart", "chromEnd", "rate")
mrmap$chr.id <- 0
for(i in 1:22) { mrmap$chr.id[which(mrmap$chrom == paste("chr", i, sep = ""))] <- i }
mrmap <- dplyr::arrange(.data = mrmap, chr.id, chromStart)
mrmap <- mrmap[-which(names(mrmap) == "chr.id")]
mrmap$rate <- mrmap$rate * 1e-8 # converts from cM/Mb to r

#######################
#
# Organizes exonic / non-exonic data
#
#######################

# to scale recombination rates within windows (Peter Ralph's formula presented in the SLiM manual)
scale_rec_rate <- function(r, n) { 
  val = (0.5 * (1 - (1 - 2 * r) ^ n)) # has underflow problems for very small r, giving 1 - 2 * r = 1
  # we adjust
  if(val == 0) {
    return(r * n)
  }
  else {
    return(val)
  }
}

human_exome_df <- as.data.frame(human_exome_gr)

# start rec. rate to zero 
human_exome_df$female_rate <- 0
human_exome_df$male_rate <- 0

# complement of the exome
non_exon_gr <- GenomicRanges::gaps(human_exome_gr)
non_exon_df <-  as.data.frame(non_exon_gr)

# to match the columns in the exome data frame
non_exon_df$ensembl_gene_id <- NA
non_exon_df$ensembl_exon_id <- NA
non_exon_df$female_rate <- 0
non_exon_df$male_rate <- 0

# merges exonic and non-exonic regions
tmp <- rbind.data.frame(human_exome_df, non_exon_df)
human_genome_df <- dplyr::arrange(.data = tmp, seqnames, start)
human_genome_gr <- makeGRangesFromDataFrame(human_genome_df, keep.extra.columns = T)

human_genome_df <- human_genome_df[,-(4:5)] # cleaning

# rec. maps Genomic Ranges
frmap_gr <- makeGRangesFromDataFrame(frmap, keep.extra.columns = T)
frmap_df <- as.data.frame(frmap_gr)
mrmap_gr <- makeGRangesFromDataFrame(mrmap, keep.extra.columns = T)
mrmap_df <- as.data.frame(mrmap_gr)

# finds gaps in rec maps relative to exonic coordinates
chr_labels <- character(length = 22)
for(i in 1:22) { chr_labels[i] <- paste("chr", i, sep = "") }
chr_ends <- c(human_genome_df[which(human_genome_df$seqnames != dplyr::lag(human_genome_df$seqnames)) - 1, 3], human_genome_df[nrow(human_genome_df), 3])
names(chr_ends) <- chr_labels
rmap_gaps <- GenomicRanges::gaps(frmap_gr, end = chr_ends)
rmap_gaps_df <- as.data.frame(rmap_gaps)
rmap_gaps_df <- rmap_gaps_df[which(rmap_gaps_df$strand == "*"),]

frmap_df <- full_join(frmap_df, rmap_gaps_df)
frmap_df <- dplyr::arrange(.data = frmap_df, seqnames, start)
frmap_df[which(is.na(frmap_df$rate)), 6] <- 0

mrmap_df <- full_join(mrmap_df, rmap_gaps_df)
mrmap_df <- dplyr::arrange(.data = mrmap_df, seqnames, start)
mrmap_df[which(is.na(mrmap_df$rate)), 6] <- 0

# female & male rec. maps share ranges...
findOverlaps(query = mrmap_gr, subject = frmap_gr, type = "equal") 

# ...so we only need to find overlaps with one
rec_ranges <- makeGRangesFromDataFrame(frmap_df, keep.extra.columns = T)

hits <- findOverlaps(query = human_genome_gr, subject = rec_ranges) 
hit_list <- as.data.frame(hits)
names(hit_list) <- c("genomic_regions", "rec_bins")

# list of genomic_regions that fall within a single recombination rate window
lookup_tbl_single <- as.data.frame(hit_list %>% group_by(genomic_regions) %>% filter(n() == 1))

# list of genomic_regions that fall within multiple recombination rate windows
lookup_tbl_multiple <- as.data.frame(hit_list %>% group_by(genomic_regions) %>% filter(n() > 1))

# first and last genomic windows covered in recombination maps
genomic_limits_in_rec_map <- c(min(c(lookup_tbl_single$genomic_regions, lookup_tbl_multiple$genomic_regions)),
                               max(c(lookup_tbl_single$genomic_regions, lookup_tbl_multiple$genomic_regions)))

frCol <- which(names(human_genome_df) == "female_rate")
mrCol <- which(names(human_genome_df) == "male_rate")

pb <- txtProgressBar(min = 0, max = nrow(human_genome_df), style = 3)
for(focal_window in 1:nrow(human_genome_df)) {
  
  setTxtProgressBar(pb, focal_window)
  
  if(focal_window >= genomic_limits_in_rec_map[1] &&
     focal_window <= genomic_limits_in_rec_map[2]) 
  { # if the genomic window matches the rec. map
    
    # if the genomic_region falls within a single recombination rate window
    if(focal_window %in% lookup_tbl_single[,1]) {
      
      row <- which(lookup_tbl_single$genomic_regions == focal_window)
      
      human_genome_df[focal_window, frCol] <- frmap_df[lookup_tbl_single[row, 2], 6]
      human_genome_df[focal_window, mrCol] <- mrmap_df[lookup_tbl_single[row, 2], 6]
    }
    
    # if it falls within multiple recombination rate windows
    else {
      
      rec_bins <- lookup_tbl_multiple[lookup_tbl_multiple$genomic_regions == focal_window, 2]
      
      human_genome_df[focal_window, frCol] <- weighted.mean(frmap_df[rec_bins, 6], frmap_df[rec_bins, 4], na.rm = T)
      human_genome_df[focal_window, mrCol] <- weighted.mean(mrmap_df[rec_bins, 6], mrmap_df[rec_bins, 4], na.rm = T)
    }
    
    # integrates rates outside exons since they will be treated as 1-bp
    if(focal_window %% 2 == 1) {
      bp <- human_genome_df[focal_window, 3] - human_genome_df[focal_window, 2] # focal window size
      human_genome_df[focal_window, frCol] <- scale_rec_rate(human_genome_df[focal_window, frCol], bp)
      human_genome_df[focal_window, mrCol] <- scale_rec_rate(human_genome_df[focal_window, mrCol], bp)
    }
  }
} 
close(pb)

# adds free rec. at chr boundaries (because SLiM technically considers the whole sequence as one chr)
human_genome_df[which(human_genome_df$seqnames != dplyr::lag(human_genome_df$seqnames)), frCol] <- 0.5
human_genome_df[which(human_genome_df$seqnames != dplyr::lag(human_genome_df$seqnames)), mrCol] <- 0.5

# to add cummulative exonic-coordinates with "1-bp" boundaries in-between
human_genome_df$exome_start <- 0
human_genome_df$exome_end <- 0
human_genome_df$exome_end[1] <- 1

exome_start_col <- which(names(human_genome_df) == "exome_start")
exome_end_col <- which(names(human_genome_df) == "exome_end")

# starts accumulating coordinates
pb <- txtProgressBar(min = 0, max = nrow(human_genome_df), style = 3)
for(i in 2:nrow(human_genome_df)) {
  
  setTxtProgressBar(pb, i)
  
  if(i %% 2 == 0) { # exon
    human_genome_df[i, exome_start_col] <- human_genome_df[i - 1, exome_end_col] 
    human_genome_df[i, exome_end_col] <- human_genome_df[i, exome_start_col] +
                                         human_genome_df[i, 3] - human_genome_df[i, 2] 
  }
  
  else { # non-coding
    human_genome_df[i, exome_start_col] <- human_genome_df[i - 1, exome_end_col] 
    human_genome_df[i, exome_end_col] <- human_genome_df[i, exome_start_col] + 1  
  }
}
close(pb)

write.table(human_genome_df, "~/Data/TIDES/human_genome/human_genome.txt", col.names = T, row.names = F, quote = F, sep = "\t")

# creates a table that will be used to format the chr column of SLiM-output VCF files 
chr.list <- character(length = 22)
for(i in 1:22) {
  chr.list[i] <- paste("chr", as.character(i), sep = "")    
}

ends.list <- numeric(length = 22)
for(i in 1:22) {
  ends.list[i] <- max(human_genome_df[which(human_genome_df$seqnames == chr.list[i]),]$exome_end)
}

chr.limits <- cbind.data.frame(chr.list, ends.list) # will be used by format_vcf.R || filter_vcfs.sh
write.table(chr.limits, "~/Data/TIDES/human_genome/chr_limits.txt", row.names = F, col.names = T, quote = F, sep = "\t")

# the condensed rec. maps (used for ABC inference)
tmp <- human_genome_df[-which(human_genome_df$female_rate == 0.5),]
write.table(tmp[c(1, exome_start_col, exome_end_col, frCol)], "~/Data/TIDES/human_genome/frmap_exome.bedgraph", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(tmp[c(1, exome_start_col, exome_end_col, mrCol)], "~/Data/TIDES/human_genome/mrmap_exome.bedgraph", col.names = F, row.names = F, quote = F, sep = "\t")

