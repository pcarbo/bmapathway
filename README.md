bmapathway
==========

This repository contains [MATLAB](www.mathworks.com/products/matlab)
source code and scripts for integrated analysis of genetic variants
and pathways in genome-wide data sets for seven complex diseases:
bipolar disorder (BD), coronary artery disease (CAD), Crohn's disease
(CD), hypertension (HT), rheumatoid arthritis (RA), type 1 diabetes
(T1D) and type 2 diabetes (T2D). These data sets are from the Wellcome
Trust Case-Control Consortium (WTCCC) studies (the initial results of
these studies were
[published in Nature in 2007](http://dx.doi.org/10.1038/nature05911)).

Running the MATLAB scripts in the [MATLAB/analysis](MATLAB/analysis)
folder should reproduce the results of our *PLoS Genetics* paper,
**"Integrated enrichment analysis of variants and pathways in
genome-wide association studies indicates central role for IL-2
signaling genes in type 1 diabetes, and cytokine signaling genes in
Crohn's disease."** For more details on the methods used, please
consult the *PLoS Genetics* paper.

This repository also contains MATLAB code implementing statistical
procedures to (1) interrogate support for enrichment of disease
associations in genome-wide data; and (2) map genetic variants
associated with disease risk. The mapping **prioritizes genetic
variants assigned to enriched gene sets, in an attempt to enhance
discovery of genes underlying complex diseases.** Our statistical
procedures are based on fitting multi-marker models of disease to the
data. We use Bayesian model averaging (BMA) in large-scale
multivariate regression to quantify support for enrichment models, and
to infer disease associations conditioned on these models.

###License

Copyright (c) 2013, Peter Carbonetto.

The bmapathway project repository by
[Peter Carbonetto](http://github.com/pcarbo) is free software: you can
redistribute it and/or modify it under the terms of the
[GNU General Public License](http://www.gnu.org/licenses/gpl.html)
as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version.

The bmapathway project repository is distributed in the hope that it
will be useful, but **without any warranty**; without even the implied
warranty of **merchantability** or **fitness for a particular
purpose**. See [LICENSE](LICENSE) for more details.

###Note about the data

The data used in the analyses are stored in the [data](data) folder:
[pathway.mat](data/pathway.mat) stores gene sets retrieved from online
pathway databases, such as [KEGG](http://www.genome.jp/kegg) and
[Reactome](http://www.reactome.org); [gene.mat](data/gene.mat) stores
information about how genes are annotated to the human genome (note
that we use version 17, or NCBI Build 35, of the Human Genome Assembly
because the data from the WTCCC disease studies are also based on this
assembly). However, we cannot make the full genotype data available
due to privacy considerations; even if we were allowed to release
these data, space restrictions on github would prevent us from storing
these files in the repository. Instead, we provide "representative"
files containing information about the genetic markers (these markers
are single nucleotide polymorphisms, or SNPs), except that the n x p
genotype matrix (where n is the number of samples, and p is the number
of SNPs) is replaced by an n x p sparse matrix, in which only a few
columns of this matrix have nonzero entries. The nonzero entries are
minor allele counts at the SNPs, in which the rows have been permuted
to preserve privacy. Thus, only the minor allele frequency of each of
these SNPs is preserved in the data we have made available. For
example, in the Crohn's disease data set, [cd.mat](data/cd.mat), the
4686 x 442,001 matrix of genotypes is replaced by a sparse matrix of
the same size, in which we have provided permuted genotypes for only
10,000 of the 442,001 SNPs.

###Overview of the MATLAB code

The [MATLAB](MATLAB) folder is organized into several subfolders. All
the MATLAB code (.m files) is found within these subfolders. There are
many files in these subfolders defining various MATLAB functions used
for our statistical analysis procedures. Here we point out the most
important folders and files, and explain when they might be useful.

+ The **[analysis](MATLAB/analysis)** folder contains the main scripts
  that run all steps of the integrated analysis for the seven
  diseases, as well as a few other functions for loading structures
  used in the analysis. All these scripts have several stages to the
  analysis. To complete the analysis, you will need to generate the
  results of each of these stages, in order. For example, the analysis
  of the Crohn's disease data set takes 11 separate steps. This
  includes computation of posterior quantities from the multi-marker
  model without pathways (Stages A and B), computation of Bayes
  factors for candidate gene sets retrieved from online pathway
  databases (Stages C and D), computation of Bayes factors for
  combinations of enriched pathways (Stages E through J), and finally
  computation of some posterior quantities conditioned on enrichment
  models with the largest Bayes factors (Stage K). All of these steps
  are implemented in the MATLAB script [cdpath.m](MATLAB/analysis/cdpath.m).

+ The **[results](MATLAB/results)** folder contains several functions
  and a script, [compileresults.m](MATLAB/results/compileresults.m),
  that compiles results from the analysis for all seven diseases, and
  generates tables and graphs for the *PLoS Genetics* paper.

+ MATLAB functions implementing our main statistical procedures are
  stored in the **[multisnp](MATLAB/multisnp)** folder. Function
  **multisnpbinhyper** runs the full variational inference procedure
  for Bayesian variable selection in logistic regression. It fits the
  multi-marker disease model to the data under the null hypothesis
  that no pathways are enriched for disease associations. Function
  **bayesfactorbin** computes the Bayes factor for a specified pathway
  annotation. It does so by fitting the multi-marker disease model to
  the data under the hypothesis that markers assigned to the pathway
  are enriched for diseases associations. Function **varpathbin**
  computes Bayes factors for a list of candidate pathways. We have
  implemented variants to each of these procedures,
  **multisnpbinhyper2**, **bayesfactorbin2** and
  **varpathbin2**. These variants are used for the modified analysis
  of the rheumatoid arthritis and type 1 diabetes data sets. These
  modifications are needed to account for the large contributions of
  MHC alleles to disease risk.

+ The **[data](MATLAB/data)** folder contains several functions and a
  script, [getwtcccdata.m](MATLAB/data/getwtcccdata.m), for acquiring
  and processing the genotype data, and for storing it in a convenient
  format for subsequent analysis steps. The genotype data was
  originally stored in files for use by the program
  [BIMBAM](http://www.bcm.edu/cnrc/mcmcmc/index.cfm?pmid=18981), but
  these files are not available here because they are large, and would
  infringe on privacy needs. Therefore, these functions are unlikely
  to be useful unless you have access to these files, or to the files
  from the original WTCCC analysis.

+ Inside the **[genes](MATLAB/genes)** folder are functions and a
  script, [getgenedata.m](MATLAB/genes/getgenedata.m), for reading
  gene information files into MATLAB. The MATLAB data file
  [gene_35.1.mat](data/gene_35.1.mat) generated by this script,
  containing information about genes corresponding to version 17 of
  the Human Genome Assembly 17 (NCBI Build 35.1), is provided in this
  repository.

+ Likewise, the **[pathways](MATLAB/pathways)** folder contains a
  number of functions and a script,
  [getpathwaydata.m](MATLAB/pathways), for reading the pathway
  annotations from the various pathway databases (Reactome, KEGG,
  PANTHER, *etc*) into MATLAB. The final MATLAB data file
  [pathway.mat](data/pathway.mat) containing the gene set annotations
  for all 3160 candidate pathways included in our analyses, is
  provided in this repository.

+ The **[mhc](MATLAB/mhc)** folder contains a couple functions used to
  specify the annotations corresponding to the major
  histocompatibility complex (MHC) and the "extended" MHC.

+ The **[stats](MATLAB/stats)** folder contains a few functions that
  calculate various statistical quantities that are useful for the
  analysis.

+ Finally, the **[misc](MATLAB/misc)** folder contains several
  miscellaneous functions that are used by the other functions.

###Downloading and using the MATLAB code

Explain how to download the code, either by downloading and unpacking
the ZIP file, or cloning the repository from the command line.

Add MATLAB directory and all subdirectories to your MATLAB path using
the addpath function in MATLAB.

You will need to download the functions implemented in the varbvs
repository, and add these functions to your path.

Code was run in MATLAB 7.10 (R2010a). Has partially been run
successfully in MATLAB 7.14 (R2012a). Has not been completely tested
in more up-to-date versions of MATLAB. Explain that 

###Credits

All MATLAB source code contained in this repository was developed by:<br>
[Peter Carbonetto](http://www.cs.ubc.ca/spider/pcarbo)<br>
Dept. of Human Genetics<br>
University of Chicago<br>
August 2013
