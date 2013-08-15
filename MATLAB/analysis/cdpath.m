% This script gives all the step I take in integrated analysis of variants
% and pathways in the WTCCC study for bipolar disorder (BD).
clear

% SCRIPT PARAMETERS.
stage = 'A';    % Which stage to run.
seed  = 7;      % Random number generator seed.
sd    = 0.1;    % Prior standard deviation of the additive effects.
da    = 100e3;  % Distance for assigning SNPs to genes (100 kb).
na    = 1e6;    % Maximum number of SNP-gene annotations.

% Initialize the random number generator. 
initrng(seed);

% LOAD GENOTYPE AND PHENOTYPE DATA.
fprintf('Loading genotype and phenotype data.\n');
load('cd.mat');

% LOAD GENE DATA.
% Load the gene data for human genome assembly 17 (NCBI Build 35) because
% the SNP data from the WTCCC case-control disease studies are based on this
% assembly.
fprintf('Loading gene data.\n');
load('gene_35.1.mat');

% LOAD PATHWAY DATA.
fprintf('Loading pathway data.\n');
load('pathway.mat');

% ASSIGN SNPS TO GENES, THEN GENES TO PATHWAYS.
fprintf('Assigning SNPs to genes and pathways.\n');
Agene  = snps2genes(chr,pos,gene,da,na);
A      = spones(Agene * pathway.genes);

% Modify the major histocompatibiltiy complex (MHC) and the annotations for
% the "extended" MHC so that they include all SNPs in each region, not just
% the ones near (known) expressed genes.
A = modifymhc(pathway,A);

fprintf('INTEGRATED ANALYSIS OF SNPs AND PATHWAYS IN CROHN''S DISEASE.\n');
fprintf('Data from study are genotypes at %d SNPs ',length(labels));
fprintf('sampled for %d cases and\n',sum(y == 1));
fprintf('%d controls.\n',sum(y == 0));
switch stage
 case 'A'

  % The first step is to conduct a multi-marker analysis of SNPs *without* any
  % information about pathways. This step is necessary for two reasons: (1)
  % to determine what regions of the genome are relevant to disease solely
  % based on genome-wide data from the case-control study so that we can
  % later appreciate what is gained by incorporating pathways into the
  % analysis; (2) to calculate posterior quantities that will later be
  % reused in computing Bayes factors for enrichment of gene sets. I begin
  % with a coarse grid for settings of the genome-wide log-odds (theta0).
  theta0 = (-6:0.25:-3)';

  % RUN JOINT ANALYSIS OF SNPs WITHOUT PATHWAYS.
  % Note that the hyperparameters used in MULTISNPBINHYPER are THETA0 and H,
  % the latter which is prior estimate of the proportion of variance
  % explained. So I need to start by converting (THETA0,SD) to H. Also note
  % that the prior on THETA0 is uniform over the range of candidate values.
  fprintf('STAGE A: joint analysis of SNPs without pathways, ');
  fprintf('given additive effects\n');
  fprintf('distribution with mean 0 and standard deviation %0.2f. ',sd);
  fprintf('Computing importance\n');
  fprintf('weights for %d candidate values of theta0, ',numel(theta0));
  fprintf('ranging from %0.2f to %0.2f.\n',min(theta0),max(theta0));
  sx = sum(var1(X));
  h  = pve(sx,sd^2,theta0);
  [logw alpha mu s eta] = multisnpbinhyper(X,y,h,theta0);
  
  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('pathway-cd-A.mat','seed','sd','theta0','logw',...
       'alpha','mu','s','eta','-v7.3');

 case 'B'

  % In this step, I continue with the multi-marker analysis of SNPs without
  % pathways. This time, I use a more fine-grained numerical approximation
  % for the average, or integral, over the genome-wide log-odds (theta0). We
  % need to restrict the range of settings for THETA0 so that this step does
  % not take too long, but we also need to choose a wide enough range of
  % settings so that we account for most of the posterior mass over the
  % genome-wide log-odds in later stages of the analysis. Based on the
  % results of Stage A and later stages examining enrichment of pathways, it
  % appears that most of the posterior mass for the genome-wide log-odds
  % lies between -6 and -3. Therefore, we can reasonably ignore settings
  % of the genome-wide log-odds outside that range.
  theta0 = (-6:0.1:-3)';

  % To speed up convergence of the coordinate ascent iterations, I initialize
  % the variational estimates of the posterior statistics using the values
  % obtained from Stage A. Specifically, I initialize the variational
  % parameters to the setting of the genome-wide log-odds with the largest
  % importance weight.
  a       = load('pathway-cd-A.mat');
  [ans i] = max(a.logw);
  alpha0  = a.alpha(:,i);
  mu0     = a.mu(:,i);
  eta0    = a.eta(:,i);

  % RUN JOINT ANALYSIS OF SNPs WITHOUT PATHWAYS.
  fprintf('STAGE B: joint analysis of SNPs without pathways, ');
  fprintf('given additive effects\n');
  fprintf('distribution with mean 0 and standard deviation %0.2f. ',sd);
  fprintf('Computing importance\n');
  fprintf('weights for %d candidate values of theta0, ',numel(theta0));
  fprintf('ranging from %0.2f to %0.2f.\n',min(theta0),max(theta0));
  sx = sum(var1(X));
  h  = pve(sx,sd^2,theta0);
  [logw alpha mu s eta] = multisnpbinhyper(X,y,h,theta0,alpha0,mu0,eta0);
  
  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('pathway-cd-B.mat','seed','sd','theta0','logw',...
       'alpha','mu','s','eta','-v7.3');

 case 'C'

  % The goal of this step of the analysis is to narrow the search for enriched
  % gene sets by obtaining rough estimates of the Bayes factors for all
  % candidate gene sets. Here, I compute the Bayes factors using a coarse
  % grid for the hyperparameters, and over a restricted range of the most
  % plausible hyperparameter settings. Later, I will use a more fine-grained
  % grid over a wider range of hyperparameters settings to refine my
  % numerical estimates of the Bayes factors for the most promising
  % candidate pathways.
  theta0 = (-6:0.5:-3.5)';
  theta  = (0:0.5:3)';

  % SPECIFY PATHWAY HYPOTHESES.
  H = getinitpaths(gene,pathway);

  % Load the posterior statistics under the null hypothesis (Stage A).
  null = getnullstats('pathway-cd-A.mat',theta0);

  % COMPUTE BAYES FACTORS FOR ENRICHMENT OF GENE SETS.
  [m n] = size(H);
  fprintf('STAGE C: integrated analysis of SNPs and pathways, ');
  fprintf('given additive effects\n');
  fprintf('distribution with mean 0 and standard deviation %0.2f, ',sd);
  fprintf('and with SNP-gene\n');
  fprintf('annotation distance of %d kb. ',da/1e3);
  fprintf('Computing Bayes factors for %d enrichment\n',n);
  fprintf('hypotheses, with %d settings for theta0 ranging ',numel(theta0));
  fprintf('from %0.2f to %0.2f, and\n',min(theta0),max(theta0));
  fprintf('with %d settings for theta ranging from %0.2f to %0.2f.\n',...
          numel(theta),min(theta),max(theta));
  [BF logw1] = varpathbin(X,H,A,y,sd^2,theta0,theta,null.logw,...
                          null.alpha,null.mu,null.eta,gene,pathway);

  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('pathway-cd-C.mat','seed','sd','da','A','H',...
       'theta0','theta','BF','logw1','-v7.3');

 case 'D'

  % Now I have identified the most promising candidates in Stage D, I compute
  % more accurate estimates of the Bayes factors for selected gene sets
  % using a finer-grained grid, and using a wider range for the
  % hyperparameters.
  theta0 = (-5:0.1:-3)';
  theta  = (0:0.1:3)';

  % SPECIFY PATHWAY HYPOTHESES.
  % I select the pathways with Bayes factors greater than 100 based on the
  % initial calculations in Stage C, and the Bayes factors for enrichment
  % of the MHC and xMHC.
  paths = [ 447 960 972 1273 1322 1391 1418 1421 1426 1700 ...
            2069 2292 2295 2423 2426 2427 2429 2591 2827 3322 3323 ];
  H     = genhmatrix(num2cell(paths),length(pathway.label));

  % Load the posterior statistics under the null hypothesis (Stage B).
  null = getnullstats('pathway-cd-B.mat',theta0);

  % COMPUTE BAYES FACTORS FOR ENRICHMENT OF GENE SETS.
  [m n] = size(H);
  fprintf('STAGE D: integrated analysis of SNPs and pathways, ');
  fprintf('given additive effects\n');
  fprintf('distribution with mean 0 and standard deviation %0.2f, ',sd);
  fprintf('and with SNP-gene\n');
  fprintf('annotation distance of %d kb. ',da/1e3);
  fprintf('Computing Bayes factors for %d enrichment\n',n);
  fprintf('hypotheses, with %d settings for theta0 ranging ',numel(theta0));
  fprintf('from %0.2f to %0.2f, and\n',min(theta0),max(theta0));
  fprintf('with %d settings for theta ranging from %0.2f to %0.2f.\n',...
          numel(theta),min(theta),max(theta));
  [BF logw1] = varpathbin(X,H,A,y,sd^2,theta0,theta,null.logw,...
                          null.alpha,null.mu,null.eta,gene,pathway);

  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('pathway-cd-D.mat','seed','sd','da','A','H',...
       'theta0','theta','BF','logw1','-v7.3');

 case 'E'

  % Now that I've identified "cytokine signaling in immune system" as the
  % pathway (or rather, this is a group of related pathways) with the
  % greatest support for enrichment of Crohn's disease associations, the
  % next step is to search for enriched pathways in combination with
  % cytokine signaling genes. I begin by computing rough estimates of the
  % Bayes factors for candidate pairs.
  theta0 = (-6:0.5:-3.5)';
  theta  = (0:0.5:3)';

  % SPECIFY PATHWAY HYPOTHESES.
  % I compute Bayes factors for the top pathway in the initial ranking from
  % Stages C and D, "cytokine signaling in immune system", combined with the
  % remaining pathways with Bayes factors greater than 10 according to the
  % initial ranking from Stage C.
  paths = [ 101 447 684 960 961 963 964 965 968 969 972 973 984 985 1012 ...
            1162 1260 1273 1289 1291 1322 1355 1357 1360 1368 1386 1391 ...
            1402 1415 1418 1421 1426 1427 1438 1448 1589 1611 1635 1700 ...
            1727 1791 1877 1898 2075 2076 2086 2093 2147 2163 2184 2222 ...
            2242 2250 2292 2295 2296 2301 2303 2304 2309 2422 2423 2426 ...
            2427 2428 2429 2525 2544 2564 2591 2653 2673 2710 2807 2811 ...
            2819 2823 2827 2941 3135 3146 3165 3192 3218 3287 ];
  H = genhmatrix(num2cell(paths),length(pathway.label));
  H(2069,:) = 1;
  
  % Load the posterior statistics under the null hypothesis (Stage A).
  null = getnullstats('pathway-cd-A.mat',theta0);

  % COMPUTE BAYES FACTORS FOR ENRICHMENT OF GENE SETS.
  [m n] = size(H);
  fprintf('STAGE E: integrated analysis of SNPs and pathways, ');
  fprintf('given additive effects\n');
  fprintf('distribution with mean 0 and standard deviation %0.2f, ',sd);
  fprintf('and with SNP-gene\n');
  fprintf('annotation distance of %d kb. ',da/1e3);
  fprintf('Computing Bayes factors for %d enrichment\n',n);
  fprintf('hypotheses, with %d settings for theta0 ranging ',numel(theta0));
  fprintf('from %0.2f to %0.2f, and\n',min(theta0),max(theta0));
  fprintf('with %d settings for theta ranging from %0.2f to %0.2f.\n',...
          numel(theta),min(theta),max(theta));
  [BF logw1] = varpathbin(X,H,A,y,sd^2,theta0,theta,null.logw,...
                          null.alpha,null.mu,null.eta,gene,pathway);

  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('pathway-cd-E.mat','seed','sd','da','A','H',...
       'theta0','theta','BF','logw1','-v7.3');

 case 'F'

  % In this step, I compute more accurate estimates of Bayes factors for
  % some of the enrichment hypotheses in which pairs of pathways are
  % enriched, using a more fine-grained grid for the hyperparameters, with
  % the range of hyperparameters selected based on the results of Stage E.
  theta0 = (-6:0.2:-3.4)';
  theta  = (1:0.25:4)';

  % SPECIFY PATHWAY HYPOTHESES.
  % I compute Bayes factors for enrichment hypotheses with Bayes factors
  % greater than 1e7 in Stage E.
  paths     = [ 684 1012 1162 1273 1421 1426 1589 1611 1635 1700 1727 ...
                1791 1877 2222 2591 2653 3135 3192 3218 3287 ];
  H         = genhmatrix(num2cell(paths),length(pathway.label));
  H(2069,:) = 1;
  
  % Load the posterior statistics under the null hypothesis (Stage B).
  null = getnullstats('pathway-cd-B.mat',theta0);

  % COMPUTE BAYES FACTORS FOR ENRICHMENT OF GENE SETS.
  [m n] = size(H);
  fprintf('STAGE F: integrated analysis of SNPs and pathways, ');
  fprintf('given additive effects\n');
  fprintf('distribution with mean 0 and standard deviation %0.2f, ',sd);
  fprintf('and with SNP-gene\n');
  fprintf('annotation distance of %d kb. ',da/1e3);
  fprintf('Computing Bayes factors for %d enrichment\n',n);
  fprintf('hypotheses, with %d settings for theta0 ranging ',numel(theta0));
  fprintf('from %0.2f to %0.2f, and\n',min(theta0),max(theta0));
  fprintf('with %d settings for theta ranging from %0.2f to %0.2f.\n',...
          numel(theta),min(theta),max(theta));
  [BF logw1] = varpathbin(X,H,A,y,sd^2,theta0,theta,null.logw,...
                          null.alpha,null.mu,null.eta,gene,pathway);

  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('pathway-cd-F.mat','seed','sd','da','A','H',...
       'theta0','theta','BF','logw1','-v7.3');

 case 'G'

  % Here I further refine the numerical estimates of the Bayes factors for
  % selected pairs of enriched pathways, and I output the posterior
  % statistics about individual SNPs conditioned on enrichment of the gene
  % sets. The range for the hyperparameters is chosen based on the results
  % of Stage F.
  theta0 = (-5.6:0.1:-3.6)';
  theta  = (1:0.1:3.5)';

  % SPECIFY PATHWAY HYPOTHESES.
  paths     = [ 684 1012 1162 1273 1421 1426 1589 1611 1635 1700 1727 ...
                1791 1877 2222 2591 2653 3135 3192 3218 3287 ];
  H         = genhmatrix(num2cell(paths),length(pathway.label));
  H(2069,:) = 1;
  
  % Load the posterior statistics under the null hypothesis (Stage B).
  null = getnullstats('pathway-cd-B.mat',theta0);

  % COMPUTE BAYES FACTORS FOR ENRICHMENT OF GENE SETS.
  [m n] = size(H);
  fprintf('STAGE G: integrated analysis of SNPs and pathways, ');
  fprintf('given additive effects\n');
  fprintf('distribution with mean 0 and standard deviation %0.2f, ',sd);
  fprintf('and with SNP-gene\n');
  fprintf('annotation distance of %d kb. ',da/1e3);
  fprintf('Computing Bayes factors for %d enrichment\n',n);
  fprintf('hypotheses, with %d settings for theta0 ranging ',numel(theta0));
  fprintf('from %0.2f to %0.2f, and\n',min(theta0),max(theta0));
  fprintf('with %d settings for theta ranging from %0.2f to %0.2f.\n',...
          numel(theta),min(theta),max(theta));
  [BF logw1 alpha mu s] = varpathbin(X,H,A,y,sd^2,theta0,theta,null.logw,...
                                     null.alpha,null.mu,null.eta,...
                                     gene,pathway);

  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('pathway-cd-G.mat','seed','sd','da','A','H',...
       'theta0','theta','BF','logw1','alpha','mu','s','-v7.3');

 case 'H'

  % In this step, I investigate support for combinations of three (3) enriched
  % pathways, again following a "greedy" approach to select combinations of
  % pathways. I take the top pair of enriched pathways, "cytokine signaling
  % in immune system" and "IL23-mediated signaling events", and compute
  % Bayes factors for these two pathways combined with the other pathways
  % with Bayes factors greater than 10 according to the initial ranking
  % (Stage C). Here I compute rough estimates for the Bayes factors using
  % a course grid for the hyperparameters.
  theta0 = (-6:0.5:-3.5)';
  theta  = (0:0.5:3)';
  
  % SPECIFY PATHWAY HYPOTHESES.
  paths = [ 101 447 684 960 961 963 964 965 968 969 972 973 984 985 1012 ...
            1162 1260 1273 1289 1291 1322 1355 1357 1360 1368 1386 1391 ...
            1402 1415 1418 1421 1426 1427 1438 1448 1589 1611 1635 1700 ...
            1727 1791 1877 1898 2075 2076 2086 2093 2147 2163 2184 2222 ...
            2242 2250 2292 2295 2296 2301 2303 2304 2309 2422 2423 2426 ...
            2427 2428 2429 2525 2544 2564 2653 2673 2710 2807 2811 2819 ...
            2823 2827 2941 3135 3146 3165 3192 3218 3287 ];
  H = genhmatrix(num2cell(paths),length(pathway.label));
  H(2069,:) = 1;
  H(2591,:) = 1;
  
  % Load the posterior statistics under the null hypothesis (Stage A).
  null = getnullstats('pathway-cd-A.mat',theta0);

  % COMPUTE BAYES FACTORS FOR ENRICHMENT OF GENE SETS.
  [m n] = size(H);
  fprintf('STAGE H: integrated analysis of SNPs and pathways, ');
  fprintf('given additive effects\n');
  fprintf('distribution with mean 0 and standard deviation %0.2f, ',sd);
  fprintf('and with SNP-gene\n');
  fprintf('annotation distance of %d kb. ',da/1e3);
  fprintf('Computing Bayes factors for %d enrichment\n',n);
  fprintf('hypotheses, with %d settings for theta0 ranging ',numel(theta0));
  fprintf('from %0.2f to %0.2f, and\n',min(theta0),max(theta0));
  fprintf('with %d settings for theta ranging from %0.2f to %0.2f.\n',...
          numel(theta),min(theta),max(theta));
  [BF logw1] = varpathbin(X,H,A,y,sd^2,theta0,theta,null.logw,...
                          null.alpha,null.mu,null.eta,gene,pathway);

  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('pathway-cd-H.mat','seed','sd','da','A','H',...
       'theta0','theta','BF','logw1','-v7.3');

 case 'I'

  % In this step, I compute more accurate estimates of Bayes factors for
  % some of the enrichment hypotheses in which three (3) pathways are
  % enriched using a more fine-grained grid for the hyperparameters, with
  % the range of hyperparameters selected based on the outcome of Stage H.
  theta0 = (-5.5:0.2:-3.5)';
  theta  = (1:0.25:4)';

  % SPECIFY PATHWAY HYPOTHESES.
  % I compute Bayes factors for enrichment hypotheses with Bayes factors
  % greater than 1e9 in Stage H.
  paths     = [ 684 1012 1162 1611 1877 2222 2429 3135 3287 ];
  H         = genhmatrix(num2cell(paths),length(pathway.label));
  H(2069,:) = 1;
  H(2591,:) = 1;
  
  % Load the posterior statistics under the null hypothesis (Stage B).
  null = getnullstats('pathway-cd-B.mat',theta0);

  % COMPUTE BAYES FACTORS FOR ENRICHMENT OF GENE SETS.
  [m n] = size(H);
  fprintf('STAGE I: integrated analysis of SNPs and pathways, ');
  fprintf('given additive effects\n');
  fprintf('distribution with mean 0 and standard deviation %0.2f, ',sd);
  fprintf('and with SNP-gene\n');
  fprintf('annotation distance of %d kb. ',da/1e3);
  fprintf('Computing Bayes factors for %d enrichment\n',n);
  fprintf('hypotheses, with %d settings for theta0 ranging ',numel(theta0));
  fprintf('from %0.2f to %0.2f, and\n',min(theta0),max(theta0));
  fprintf('with %d settings for theta ranging from %0.2f to %0.2f.\n',...
          numel(theta),min(theta),max(theta));
  [BF logw1] = varpathbin(X,H,A,y,sd^2,theta0,theta,null.logw,...
                          null.alpha,null.mu,null.eta,gene,pathway);

  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('pathway-cd-I.mat','seed','sd','da','A','H',...
       'theta0','theta','BF','logw1','-v7.3');

 case 'J'

  % In Stage J, I further refine the numerical estimates of the Bayes factors
  % for enrichment hypotheses in which three (3) pathways are enriched, and I
  % output posterior statistics about individual SNPs conditioned on
  % enrichment of these gene sets. The range of hyperparameters is based
  % on the results of Stage I.
  theta0 = (-5.3:0.1:-3.7)';
  theta  = (1.5:0.1:3.3)';

  % SPECIFY PATHWAY HYPOTHESES.
  % I compute Bayes factors for enrichment hypotheses with Bayes factors
  % greater than *1e10* in Stage H. (According to the results of Stage H,
  % the next largest Bayes factor smaller than 1e10 is 40 times smaller than
  % the largest one Bayes factor.)
  paths     = [ 684 1012 1162 1611 1877 2222 3135 3287 ];
  H         = genhmatrix(num2cell(paths),length(pathway.label));
  H(2069,:) = 1;
  H(2591,:) = 1;
  
  % Load the posterior statistics under the null hypothesis (Stage B).
  null = getnullstats('pathway-cd-B.mat',theta0);

  % COMPUTE BAYES FACTORS FOR ENRICHMENT OF GENE SETS.
  [m n] = size(H);
  fprintf('STAGE J: integrated analysis of SNPs and pathways, ');
  fprintf('given additive effects\n');
  fprintf('distribution with mean 0 and standard deviation %0.2f, ',sd);
  fprintf('and with SNP-gene\n');
  fprintf('annotation distance of %d kb. ',da/1e3);
  fprintf('Computing Bayes factors for %d enrichment\n',n);
  fprintf('hypotheses, with %d settings for theta0 ranging ',numel(theta0));
  fprintf('from %0.2f to %0.2f, and\n',min(theta0),max(theta0));
  fprintf('with %d settings for theta ranging from %0.2f to %0.2f.\n',...
          numel(theta),min(theta),max(theta));
  [BF logw1 alpha mu s] = varpathbin(X,H,A,y,sd^2,theta0,theta,null.logw,...
                                     null.alpha,null.mu,null.eta,...
                                     gene,pathway);

  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('pathway-cd-J.mat','seed','sd','da','A','H',...
       'theta0','theta','BF','logw1','alpha','mu','s','-v7.3');

 case 'K'

  % In this step, I compute posterior statistics about individual SNPs
  % conditioned on selected enrichment hypotheses in which one pathway is
  % enriched.
  theta0 = (-4.5:0.1:-3.5)';
  theta  = (1:0.1:3)';

  % SPECIFY PATHWAY HYPOTHESES.
  % I select the top pathways based on the results of Stage D.
  H = genhmatrix(num2cell([ 1421 1426 2069 ]),length(pathway.label));

  % Load the posterior statistics under the null hypothesis (Stage B).
  null = getnullstats('pathway-cd-B.mat',theta0);

  % COMPUTE BAYES FACTORS FOR ENRICHMENT OF GENE SETS.
  [m n] = size(H);
  fprintf('STAGE K: integrated analysis of SNPs and pathways, ');
  fprintf('given additive effects\n');
  fprintf('distribution with mean 0 and standard deviation %0.2f, ',sd);
  fprintf('and with SNP-gene\n');
  fprintf('annotation distance of %d kb. ',da/1e3);
  fprintf('Computing Bayes factors for %d enrichment\n',n);
  fprintf('hypotheses, with %d settings for theta0 ranging ',numel(theta0));
  fprintf('from %0.2f to %0.2f, and\n',min(theta0),max(theta0));
  fprintf('with %d settings for theta ranging from %0.2f to %0.2f.\n',...
          numel(theta),min(theta),max(theta));
  [BF logw1 alpha mu s] = varpathbin(X,H,A,y,sd^2,theta0,theta,null.logw,...
                                     null.alpha,null.mu,null.eta,...
                                     gene,pathway);

  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('pathway-cd-K.mat','seed','sd','da','A','H',...
       'theta0','theta','BF','logw1','alpha','mu','s','-v7.3');
end
