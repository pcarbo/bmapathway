% This script gives all the step I take in integrated analysis of variants
% and pathways in the WTCCC study for coronary artery disease (CAD).
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
load('cad.mat');

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

fprintf('INTEGRATED ANALYSIS OF SNPs AND PATHWAYS IN ');
fprintf('CORONARY ARTERY DISEASE.\n');
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
  save('pathway-cad-A.mat','seed','sd','theta0',...
       'logw','alpha','mu','s','eta','-v7.3');

 case 'B'

  % In this step, I continue with the multi-marker analysis of SNPs without
  % pathways. This time, I use a more fine-grained numerical approximation
  % for the average, or integral, over the genome-wide log-odds (theta0). We
  % need to restrict the range of settings for THETA0 so that this step does
  % not take too long, but we also need to choose a wide enough range of
  % settings so that we account for most of the posterior mass over the
  % genome-wide log-odds in later stages of the analysis. Based on the
  % results of Stage A, it appears that most of the posterior mass for the
  % genome-wide log-odds lies between -6 and -3.5. Therefore, we can
  % reasonably ignore settings of the genome-wide log-odds outside that
  % range.
  theta0 = (-6:0.1:-3.5)';

  % To speed up convergence of the coordinate ascent iterations, I initialize
  % the variational estimates of the posterior statistics using the values
  % obtained from Stage A. Specifically, I initialize the variational
  % parameters to the setting of the genome-wide log-odds with the largest
  % importance weight.
  a       = load('pathway-cad-A.mat');
  [ans i] = max(a.logw);
  alpha0  = a.alpha(:,i);
  mu0     = a.mu(:,i);
  eta0    = a.eta(:,i);

  % RUN JOINT ANALYSIS OF SNPs WITHOUT PATHWAYS.
  % Note that the hyperparameters used in MULTISNPBINHYPER are THETA0 and H,
  % the latter which is prior estimate of the proportion of variance
  % explained. So I need to start by converting (THETA0,SD) to H. Also note
  % that the prior on THETA0 is uniform over the range of candidate values.
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
  save('pathway-cad-B.mat','seed','sd','theta0',...
       'logw','alpha','mu','s','eta','-v7.3');

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
  null = getnullstats('pathway-cad-A.mat',theta0);

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
  save('pathway-cad-C.mat','seed','sd','da','A','H',...
       'theta0','theta','BF','logw1','-v7.3');

 case 'D'

  % Now I have identified the most promising candidates in Stage C, I compute
  % more accurate estimates of the Bayes factors for selected gene sets
  % using a finer-grained grid, and using a wider range for the
  % hyperparameters.
  theta0 = (-6:0.2:-3.6)';
  theta  = (0:0.25:5)';

  % SPECIFY PATHWAY HYPOTHESES.
  % I select the pathways with Bayes factors greater than 10 according on the
  % initial calculations in Stage C.
  paths = [  323   327   328  1243  1250  1254  1272  1277  1316 ...
	    1347  1362  1365  1369  1382  1390  1398  1405  1407 ... 
	    1412  1417  1419  1429  1508  1520  1525  1748  2053 ...
	    2118  2197  2202  2210  2225  2326  2420  2618  3084 ...
	    3085  3150  3151  3153  3154  3155  3159  3161  3205 ...
	    3214  3215 ];
  H     = genhmatrix(num2cell(paths),length(pathway.label));

  % Load the posterior statistics under the null hypothesis (Stage B).
  null = getnullstats('pathway-cad-B.mat',theta0);

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
  save('pathway-cad-D.mat','seed','sd','da','A','H',...
       'theta0','theta','BF','logw1','-v7.3');

 case 'E'

  % In this step, I further refine the numerical estimates of the Bayes
  % factors for the top-ranked gene sets.
  theta0 = (-6:0.1:-3.8)';
  theta  = (1:0.1:4.5)';

  % SPECIFY PATHWAY HYPOTHESES.
  % I select the pathways with Bayes factors greater than 50 based on the
  % results of the previous step (Stage D).
  paths = [ 323 1243 1508 1748 2118 3084 3205 ];
  H     = genhmatrix(num2cell(paths),length(pathway.label));

  % Load the posterior statistics under the null hypothesis (Stage B).
  null = getnullstats('pathway-cad-B.mat',theta0);

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
  save('pathway-cad-E.mat','seed','sd','da','A','H',...
       'theta0','theta','BF','logw1','-v7.3');
end