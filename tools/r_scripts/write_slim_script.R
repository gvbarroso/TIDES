# Date created: 09/11/2020
# Author: Gustavo V. Barroso
# This script writes SLiM scripts to generate the parents of our simulated trios.
# TIDES finishes the simulation internally, switching to viability selection

# first we create df with parameter combinations
#s = c(-1e-2, -1e-3, -1e-4)
#h = c(0, 0.5)
#replicate = c(1:10)

#library(tidyverse)
#df <- crossing(s, h, replicate)
#write.table(df, "slim_simulations/positive_sel/params_df_positive.txt", col.names = F, row.names = T, quote = F, sep = "\t")

# we will params_df.txt as argument for this script using bash:

# while IFS='\t' read -r -a vals; do Rscript --vanilla write_slim_script.R ${vals}; done < params_df.txt

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

idx <- as.numeric(args[1]) # which sim index
s <- as.numeric(args[2]) # which selection coefficient
h <- as.numeric(args[3]) # which dominance coefficient
r <- as.numeric(args[4]) # which replicate

dir.create(paste("script_", idx, sep = ""))

# table generated with get_hg_data.R
exome <- read.table("~/Data/TIDES/human_genome/human_genome.txt", header = T)

# starts writing SLiM script
sink(paste("script_", idx, "/exome_s_", as.character(s), "_h_", as.character(h), "_rep_", r,
           "_deCODE_maps_gravel_demo.slim", sep = ""))

cat("initialize() {\n") # this is a Wright-Fisher simulation

cat(paste("setSeed(", as.character(r), ");\n", sep = "")) # for reproducibility

cat(paste("defineConstant(\"h\",", as.character(h), ");\n", sep = "")) 
cat(paste("defineConstant(\"s\",", as.character(s), ");\n", sep = ""))

# advantageous mutations are 10X STRONGER and are 20X LESS frequent than deleterious mutations:
cat("initializeMutationRate(8.64e-9 * 2.31 / 3.1);\n") # estimate of the human non-synonymous (NS) mutation rate
cat("initializeMutationType(\"m1\", h, \"f\", s);\n") # deleterious mutations
cat("initializeMutationType(\"m2\", h, \"f\", -10 * s);\n") # advantageous mutations 
cat("initializeGenomicElementType(\"g1\", c(m1, m2), c(0.99, 0.01));\n") # exon

cat("initializeSex(\"A\");\n\n")

# female recombination maps
cat("ratesF = c(")
for(i in 1:nrow(exome))
{
  cat(paste("asFloat(", as.character(exome$female_rate[i]), ")", sep = ""))
  if(i < nrow(exome)) { cat(",\n") }
}
cat(");\n\n")

cat("endsF = c(")
for(i in 1:nrow(exome))
{
  cat(paste("asInteger(", as.character(exome$exome_end[i]), ")", sep = ""))
  if(i < nrow(exome)) { cat(",\n") }
}
cat(");\n\n")

cat("initializeRecombinationRate(ratesF, endsF, \"F\");\n\n")

# male recombination maps
cat("ratesM = c(")
for(i in 1:nrow(exome))
{
  cat(paste("asFloat(", as.character(exome$male_rate[i]), ")", sep = ""))
  if(i < nrow(exome)) { cat(",\n") }
}
cat(");\n\n")

cat("endsM = c(")
for(i in 1:nrow(exome))
{
  cat(paste("asInteger(", as.character(exome$exome_end[i]), ")", sep = ""))
  if(i < nrow(exome)) { cat(",\n") }
}
cat(");\n\n")

cat("initializeRecombinationRate(ratesM, endsM, \"M\");\n\n")

# exome coords
exome_start_col <- which(names(exome) == "exome_start")
exome_end_col <- which(names(exome) == "exome_end")

for(i in seq(from = 2, to = nrow(exome), by =  2))
{
  cat(paste("initializeGenomicElement(g1,asInteger(", as.character(exome[i, exome_start_col]), "),", "asInteger(", as.character(exome[i, exome_end_col]), "));\n", sep = ""))
}
cat("}\n\n") # closing the "initialize" brackets

# Demography (Gravel model without Asia and with additional generations to increase sample (pop.) size)

cat("1 { sim.addSubpop(\"p1\", 7310); }\n\n")

cat("52080 { p1.setSubpopulationSize(14474); }\n\n")

cat("55960 {\n")
cat("sim.addSubpopSplit(\"p2\", 1861, p1);\n")
cat("p1.setMigrationRates(p2, 15e-5);\n")
cat("p2.setMigrationRates(p1, 15e-5);\n")
cat("}\n\n")

cat("57080 {\n")
cat("p2.setSubpopulationSize(1032);\n")
cat("p1.setMigrationRates(p2, 2.5e-5);\n")
cat("p2.setMigrationRates(p1, 2.5e-5);\n")
cat("}\n\n")

cat("57080:58150 {\n")
cat("t = sim.generation - 57080;\n")
cat("p2_size = round(1032 * exp(0.0038 * t));\n")
cat("p2.setSubpopulationSize(asInteger(p2_size));\n")
cat("}\n\n")

# Output entire parental generation (parents will undergo viability selection inside TIDES)
cat("58151 early() {\n")

cat("potentialMothers = NULL;\n")
cat("potentialFathers = NULL;\n")

cat("for(indv in p2.individuals) {\n")
cat("if(indv.sex == \"F\") {\n")
cat("potentialMothers = c(potentialMothers, indv);\n")
cat("}\n")
cat("else {\n")
cat("potentialFathers = c(potentialFathers, indv);\n")
cat("}\n")
cat("}\n")
cat("potentialMothers.genomes.outputVCF(\"potential_mothers.vcf\");\n")
cat("potentialFathers.genomes.outputVCF(\"potential_fathers.vcf\");\n")
cat("}\n\n")

cat("58151 late() { sim.simulationFinished(); }")

sink()
