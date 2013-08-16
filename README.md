bmapathway
==========

This repository contains [MATLAB](www.mathworks.com/products/matlab)
source code and scripts implementing integrated analysis of genetic
variants and pathways in genome-wide association studies (GWAS) of
seven complex diseases: bipolar disorder (BD), coronary artery disease
(CAD), Crohn's disease (CD), hypertension (HT), rheumatoid arthritis
(RA), type 1 diabetes (T1D) and type 2 diabetes (T2D). Running these
MATLAB scripts should reproduce the results on enriched pathways and
disease associations reported in the PLoS Genetics paper, *Integrated
enrichment analysis of variants and pathways in genome-wide
association studies indicates central role for IL-2 signaling genes in
type 1 diabetes, and cytokine signaling genes in Crohn's disease.*

Also included in this repository is MATLAB code implementing
statistical procedures to (1) interrogate pathways for enrichment of
disease associations in genome-wide data; and (2) map genetic variants
associated with susceptibility with disease, and specifically it
*prioritizes variants assigned to enriched pathways in an attempt to
enhance discovery of genes underlying complex diseases.* These
procedures are based on fitting multi-marker models of disease to the
data. We use Bayesian model averaging (BMA) for large-scale
multivariate regression to quantify support for models of enriched
pathways, and to infer disease associations across the genome.

Explain that data cannot be made available to the public. Even if we
were allowed to release the data, github would not permit us to store
all the data because it would take up several Gb of space.

Code was run in MATLAB 7.10 (R2010a). Has partially been run
successfully in MATLAB 7.14 (R2012a). Has not been completely tested
in more up-to-date versions of MATLAB.

Requires varbvs code!

Add MATLAB directory and all subdirectories to your MATLAB path using
the addpath function in MATLAB.

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

###Overview of the MATLAB source code

Text goes here.

+ [filename](pathtofile) Decription of file goes here.

###Credits

All MATLAB source code contained in this repository was developed by:<br>
[Peter Carbonetto](http://www.cs.ubc.ca/spider/pcarbo)<br>
Dept. of Human Genetics<br>
University of Chicago<br>
August 2013
