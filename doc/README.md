# Running TIDES

This short manual describes TIDES from a practical standpoint. 
The 'example' directory contains step-by-step instructions for using TIDES in 'example_run.Rmd'.

To run the software, the user must fill in the options file 'opt.bpp' and execute the program in the command line:

	TIDES params=opt.bpp

## A note about computing clusters

Before we describe the available options, here is some useful information about running TIDES on a cluster (credit to Shao-Ching Huang @UCLA).
The user must execute the program with the same versions of the libraries that were used to compile it.
For example, if TIDES was compiled with gcc 9.3.0, and that is not the default gcc version on the cluster, then one must load gcc 9.3.0 before executing TIDES:
 
	module load gcc/9.3.0
	./TIDES params=opt.bpp

It may be that "module load" works on a interactive session of the cluster\ but not inside the script that is used to submit the job to computing nodes. In that case, here is the structure of a cluster script that did the trick for me:


	#!/bin/sh

	[submission parameters such as memory and number of CPU requested]

	source /u/local/Modules/default/init/bash

	module load gcc/9.3.0

	./TIDES params=opt.bpp
	

## Filling in the options file

With that out of the way, we move on to TIDES itself.
Here is how the user must specify the options.

### Input files

Input sequences must be phased and provided as VCF files for mothers, fathers and children separately. Under the current version, **these must contain only sites believed to be under natural selection**, e.g., non-synonymous derived variants within exons. The trio relationships are provided in the trio_ids_file option, where each line contains three TAB-separated elements that must be provided in order: child_id, mother_id, father_id (these are the ID's given in the header of the corresponding VCF files). Trios may be overlapping, i.e., the data set may include siblings and half-siblings.

The last piece of input files are (sex-specific) recombination maps that are used in the meiosis simulations. The maps should be in BEDGRAPH format and must present rates in units of r -- cross-over per base pair per generation -- rather than in cM/Mb. If only a gender neutral map is available, it should simply be duplicated in the options file, provided once for each sex. If a recombination map is not available at all for the species of interest, the user must still write a simple BEDGRAPH file, where each line contains the start and end of each chromosome, and the rate is an estimate of the genome-wide average r. Recombination maps should cover all sites in the VCF files and may include regions outside sites in the VCF files.

Note that the VCF and the BEDGRAPH files must match their notation, i.e., if chromosomes are labeled 'chr1', 'chr2', ..., in the VCF files, they should be labeled the same in the BEDGRAPH files. Also, if genomic coordinates restart at each chromosome in the VCF files (e.g, 1-250,000,000 for chr1, then 1-225,000,000 for chr2 and so on), they must also do so in the BEDGRAPH files. Finally, it is allowed to provide absolute genomic coordinates that do not restart at each chromosome, but in any case the largest position allowed for a site is 4,294,967,295 (which is still > 1Gb longer the the human genome).

### Other parameters

In the params file (opt.bpp) users must provide a range of values for the uniform priors of s and h used in the pilot simulations. Here, s_interval and h_interval are comma-separated real numbers within parenthesis, representing the lower and upper bound of parameter values that will be explored by TIDES. **This also makes it convenient to fit different models to the data**. For example, to fit a 'recessive' model, one can set h_interval = (0, 0). Likewise, an 'additive' model is achieved by setting h_interval = (0.5, 0.5). A neutral model, on the other hand, can be specified by s_interval = (100, 100) (since the prior on s during the pilot run is log-uniform). By fitting different models to the same data, one can perform model selection.

TIDES handles de novo mutations in one of two ways. If a mutation rate per nucleotide per generation ('mu' option) is input alongside the number of target sites in the exome ('seq_len' option), de novo mutations are included in the simulation of zygotes. Note that in this case 'mu' is supposed to be the rate of **non-synonymous** mutations. Alternatively, if both options are left blank, TIDES identifies SNPs with zero frequency in the parents and removes them from the children’s genomes, in which case they are not part of the simulation of zygotes.

Finally, it is possible to specify the number of pilot simulations (where parameters are drawn from flat priors), the number of second round simulations (where parameters are drawn from the pilot posteriors), the number of meiosis (simulated zygotes) for each trio, as well as the number of accepted pilot simulations (how many posterior samples are used to build the prior for the second round of simulations). A brief discussion about changing the default value of these options is presented in the 'General Comments' section below.

## Output files

TIDES will output two compressed text files (gzipped):

- sim_params_dist.gz presents parameter values from all simulations from the second round
- sum_stats_dist.gz presents summary statistics from all simulations from the second round

These are used as input by '/tools/r_scripts/abc_filter.Rmd' in order to perform model selection and parameter inference. That script contains further comments about the procedure.

## Performance

TIDES is fairly efficient. On a Linux machine with 16 threads and 3.5 Ghz clock speed, TIDES can simulate ~ 40,000 zygotes (with a human-like exome) per second using high-resolution recombination maps. On a data set of 50,000 family trios and 5,000,000 simulated zygotes, it will simulate 11,000,000 rounds of meiosis + viability selecton in ~ 24 hours. The program is highly parallelizable such that access to computing clusters with >> 16 threads can dramatically increase its speed.

## General Comments

TIDES provides default values for all options discussed above (see comments in 'opt.bpp'). Although an exhaustive exploration of different combinations was not carried out, the default values showed robust results in our analyses. However, we encourage users to test different values that may  best suits their needs, emphasizing the trade-off between accuracy and run-time. Ideally, running TIDES with a lower number of simulations, higher acceptance rates in the pilot simulations or less simulated zygotes per trio should be validated by benchmarking on simulated data to test statistical power for both model selection and parameter inference. 