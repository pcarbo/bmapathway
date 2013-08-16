bmapathway
==========

This repository contains [MATLAB](www.mathworks.com/products/matlab)
source code and scripts for integrated analysis of genetic variants
and pathways in genome-wide association studies (GWAS) of seven
complex diseases: bipolar disorder (BD), coronary artery disease
(CAD), Crohn's disease (CD), hypertension (HT), rheumatoid arthritis
(RA), type 1 diabetes (T1D) and type 2 diabetes (T2D). These data are
from the Wellcome Trust Case-Control Consortium (WTCCC) studies,
originally published in a
[Nature paper](http://dx.doi.org/10.1038/nature05911) in 2007. Running
these MATLAB scripts should reproduce the results given in the *PLoS
Genetics* paper, **Integrated enrichment analysis of variants and
pathways in genome-wide association studies indicates central role for
IL-2 signaling genes in type 1 diabetes, and cytokine signaling genes
in Crohn's disease.** For more details on these methods, consult the
*PLoS Genetics* paper.

This repository also contains MATLAB code implementing statistical
procedures to (1) interrogate support for enrichment of disease
associations in genome-wide data; and (2) map genetic variants
associated with disease risk, including *prioritization of variants
assigned to enriched gene sets, in an attempt to enhance discovery of
genes underlying complex diseases.* These statistical procedures are
based on fitting multi-marker models of disease to the data. We use
Bayesian model averaging (BMA) in large-scale multivariate regression
to quantify support for enrichment models, and to infer disease
associations conditioned on these models.

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

###A note about the data

Explain that data cannot be made available to the public. Even if we
were allowed to release the data, github would not permit us to store
all the data because it would take up several Gb of space.

###Overview of the MATLAB code

The [MATLAB](MATLAB) folder is organized into several subfolders, and
all the MATLAB code (in .m files) is within these subfolders. There
are a lot of files in these subfolders that define various MATLAB
functions. Here we point out the most important folders and files, and
explain when they might be useful.

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

+ 

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
