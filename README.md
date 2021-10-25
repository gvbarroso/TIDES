# TIDES
Trio-based inference of dominance and selection


TIDES is a statistical model for inferring dominance and the strength of natural selection from trio sequence data. In a nutshell, TIDES is a fast simulator of meiosis and natural selection. Model critique and parameter inference are performed via approximate Bayesian computation in a step following up the simulations, using dedicated R packages (see below).

The INSTALL.txt file provides instructions on how to compile TIDES.
The 'doc/' directory contains a (hopefully useful) README file with practical considerations, an annotated options file discribing supported options, as well as one running example.
The 'tools/r_scripts' directory includes R scripts for analyzing and plotting results, as well as for reproducing the simulation results from the original article.
This project uses C++ code modified after the Bio++ Core libraries (https://github.com/BioPP/bpp-core)

To cite TIDES, use: https://www.biorxiv.org/content/10.1101/2021.10.08.463705v1
