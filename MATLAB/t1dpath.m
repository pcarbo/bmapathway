% This script gives all the step I take in integrated analysis of variants
% and pathways in the WTCCC study for type 1 diabetes (T1D).
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
load('/tmp/pcarbo/t1d.mat');

% Sets of SNPs corresponding to the MHC region (mhcsnps), and the entire
% genome except for the MHC (othersnps).
isinmhc   = inmhcregion(chr,pos);
mhcsnps   = find(isinmhc);
othersnps = find(~isinmhc);

% LOAD GENE DATA.
% Load the gene data for human genome assembly 17 (NCBI Build 35) because
% the SNP data from the WTCCC case-control disease studies are based on this
% assembly.
fprintf('Loading gene data.\n');
load('/tmp/pcarbo/gene_35.1.mat');

% LOAD PATHWAY DATA.
fprintf('Loading pathway data.\n');
load('/tmp/pcarbo/pathway.mat');

% ASSIGN SNPS TO GENES, THEN GENES TO PATHWAYS.
fprintf('Assigning SNPs to genes and pathways.\n');
Agene  = snps2genes(chr,pos,gene,da,na);
A      = spones(Agene * pathway.genes);

% Modify the major histocompatibiltiy complex (MHC) and the annotations for
% the "extended" MHC so that they include all SNPs in each region, not just
% the ones near (known) expressed genes.
A = modifymhc(pathway,A);

fprintf('INTEGRATED ANALYSIS OF SNPs AND PATHWAYS IN TYPE 1 DIABETES.\n');
fprintf('Data from study are genotypes at %d SNPs ',length(labels));
fprintf('sampled for %d cases and\n',sum(y == 1));
fprintf('%d controls.\n',sum(y == 0));
switch stage
 case 'A'

  % The first step is to conduct a multi-marker analysis of SNPs without the
  % pathways, and by treating all SNPs in the genome in the same way. As we
  % will see, there are good reasons for treating SNPs in the MHC region
  % differently in analysis of the type 1 diabetes study, but here it
  % is useful to take this first step to estimate parameters in the
  % variational approximation that will be used in later steps of the
  % analysis. I begin with a coarse grid for settings of the genome-wide
  % log-odds parameter.
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
  save('/tmp/pcarbo/pathway-t1d-A.mat','seed','sd','theta0','logw',...
       'alpha','mu','s','eta','-v7.3');

 case 'B'

  % In this step, I redo the same analysis without pathways, except that I
  % adopt a two-stage approach in which I ignore SNPs within the major
  % histocompatibility complex (MHC) when analyzing SNPs outside the MHC.
  % Again, I use a coarse grid for settings of the genome-wide log-odds
  % parameter (theta0).
  theta0 = (-6:0.25:-3)';

  % Import the estimates of the variational parameters "eta" from Stage A
  % of the analysis.
  a   = load('/tmp/pcarbo/pathway-t1d-A.mat');
  eta = a.eta;

  % RUN JOINT ANALYSIS OF SNPs WITHOUT PATHWAYS.
  % Note that the hyperparameters used in MULTISNPBINHYPER are THETA0 and H,
  % the latter which is prior estimate of the proportion of variance
  % explained. So I need to start by converting (THETA0,SD) to H. Also note
  % that the prior on THETA0 is uniform over the range of candidate values.
  fprintf('STAGE B: joint analysis of SNPs without pathways, ');
  fprintf('given additive effects\n');
  fprintf('distribution with mean 0 and standard deviation %0.2f, ',sd);
  fprintf('and analyzing SNPs\n');
  fprintf('in the MHC (%d SNPs) ',length(mhcsnps));
  fprintf('conditioned on SNPs outside the MHC (%d SNPs).\n',...
          length(othersnps));
  fprintf('Computing importance weights for %d candidate values of ',...
          numel(theta0));
  fprintf('theta0, ranging\n');
  fprintf('from %0.2f to %0.2f.\n',min(theta0),max(theta0));
  sx = sum(var1(X));
  h  = pve(sx,sd^2,theta0);
  [logw alpha mu s] = multisnpbinhyper2(X,y,othersnps,mhcsnps,h,theta0,eta);
  
  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('/tmp/pcarbo/pathway-t1d-B.mat','seed','sd','theta0',...
       'mhcsnps','othersnps','logw','alpha','mu','s','eta',...
       '-v7.3');

 case 'C'

  % In this step, I redo the multi-marker analysis of SNPs, as in Stage A, but
  % this time using a more fine-grained numerical approximation for the
  % average, or integral, over the genome-wide log-odds (theta0). We need to
  % restrict the range of settings for THETA0 so that this step does not
  % take too long, but we also need to choose a wide enough range of
  % settings so that we account for most of the posterior mass over the
  % genome-wide log-odds in later stages of the analysis. Based on the
  % results of Stage A, and from later stages examining enrichment of
  % pathways, it appears that most of the posterior mass for the genome-wide
  % log-odds lies between -6 and -3. Therefore, we can reasonably ignore
  % settings of the genome-wide log-odds outside that range.
  theta0 = (-6:0.1:-3)';

  % To speed up convergence of the coordinate ascent iterations, I initialize
  % the variational estimates of the posterior statistics using the values
  % obtained from Stage A. Specifically, I initialize the variational
  % parameters to the setting of the genome-wide log-odds with the largest
  % importance weight.
  a       = load('/tmp/pcarbo/pathway-t1d-A.mat');
  [ans i] = max(a.logw);
  alpha0  = a.alpha(:,i);
  mu0     = a.mu(:,i);
  eta0    = a.eta(:,i);

  % RUN JOINT ANALYSIS OF SNPs WITHOUT PATHWAYS.
  % Note that the hyperparameters used in MULTISNPBINHYPER are THETA0 and H,
  % the latter which is prior estimate of the proportion of variance
  % explained. So I need to start by converting (THETA0,SD) to H. Also note
  % that the prior on THETA0 is uniform over the range of candidate values.
  fprintf('STAGE C: joint analysis of SNPs without pathways, ');
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
  save('/tmp/pcarbo/pathway-t1d-C.mat','seed','sd','theta0','logw',...
       'alpha','mu','s','eta','-v7.3');

 case 'D'

  % In this step, I analyze the SNPs just as I did in Stage B, in which we
  % ignore MHC SNPs when analyzing SNPs outside the MHC, but now using a
  % more fine-grained numerical approximation for the averages, or
  % integrals, over the genome-wide log-odds hyperparameter (theta0), like
  % I did in Stage C.
  theta0 = (-6:0.1:-3)';

  % Import the estimates of the variational parameters "eta" from Stage C
  % of the analysis.
  a   = load('/tmp/pcarbo/pathway-t1d-C.mat');
  eta = a.eta;

  % RUN JOINT ANALYSIS OF SNPs WITHOUT PATHWAYS.
  % Note that the hyperparameters used in MULTISNPBINHYPER are THETA0 and H,
  % the latter which is prior estimate of the proportion of variance
  % explained. So I need to start by converting (THETA0,SD) to H. Also note
  % that the prior on THETA0 is uniform over the range of candidate values.
  fprintf('STAGE D: joint analysis of SNPs without pathways, ');
  fprintf('given additive effects\n');
  fprintf('distribution with mean 0 and standard deviation %0.2f, ',sd);
  fprintf('and analyzing SNPs\n');
  fprintf('in the MHC (%d SNPs) ',length(mhcsnps));
  fprintf('conditioned on SNPs outside the MHC (%d SNPs).\n',...
          length(othersnps));
  fprintf('Computing importance weights for %d candidate values of ',...
          numel(theta0));
  fprintf('theta0, ranging\n');
  fprintf('from %0.2f to %0.2f.\n',min(theta0),max(theta0));
  sx = sum(var1(X));
  h  = pve(sx,sd^2,theta0);
  [logw alpha mu s] = multisnpbinhyper2(X,y,othersnps,mhcsnps,h,theta0,eta);
  
  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('/tmp/pcarbo/pathway-t1d-D.mat','seed','sd','theta0',...
       'mhcsnps','othersnps','logw','alpha','mu','s','eta',...
       '-v7.3');

 case 'E'

  % The goal of this step of the analysis is to narrow the search for enriched
  % gene sets by obtaining rough estimates of the Bayes factors for all
  % candidate gene sets. Here, I compute the Bayes factors using a coarse
  % grid for the hyperparameters, and over a restricted range of the most
  % plausible hyperparameter settings. Later, I will use a more fine-grained
  % grid over a wider range of hyperparameters settings to refine my
  % numerical estimates of the Bayes factors for the most promising
  % candidate pathways. As in Stages B and D, I ignore MHC SNPs for the
  % multi-marker analysis of SNPs outside the MHC.
  theta0 = (-6:0.5:-4)';
  theta  = (0:0.5:3)';

  % SPECIFY PATHWAY HYPOTHESES.
  H = getinitpaths(gene,pathway);

  % Load the posterior statistics under the null hypothesis (Stage B).
  null = getnullstats('/tmp/pcarbo/pathway-t1d-B.mat',theta0);

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
  [BF logw1] = varpathbin2(X,H,A,y,~isinmhc,sd^2,theta0,theta,...
                           null.logw,null.alpha,null.mu,null.eta,...
                           gene,pathway);

  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('/tmp/pcarbo/pathway-t1d-E.mat','seed','sd','da','A','H',...
       'theta0','theta','BF','logw1','-v7.3');

 case 'F'

  % In the previous step, we found by far the strongest support for enrichment
  % of the MHC and the "extended" MHC (which I will abbreviate by
  % "xMHC"). In this step, I reassess the support for enrichment of the MHC
  % and xMHC (plus a few other pathways) using a wider range of settings for
  % the enrichment parameter.
  theta0 = (-6:0.25:-3)';
  theta  = (0:0.5:5)';

  % SPECIFY PATHWAY HYPOTHESES.
  paths = [ 3322 3323 3145 3146 3195 ];
  H     = genhmatrix(num2cell(paths),length(pathway.label));

  % Load the posterior statistics under the null hypothesis (Stage B).
  null = getnullstats('/tmp/pcarbo/pathway-t1d-B.mat',theta0);

  % COMPUTE BAYES FACTORS FOR ENRICHMENT OF GENE SETS.
  % Note that there is no need to adopt the modified analysis here (as
  % implemented in function VARPATHBIN2) because I am already analyzing SNPs
  % in the pathway---that is, SNPs in the MHC region---only after first
  % estimating posterior statistics of the SNPs outside the MHC, which were
  % obtained under the null hypothesis.
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
  save('/tmp/pcarbo/pathway-t1d-F.mat','seed','sd','da','A','H',...
       'theta0','theta','BF','logw1','-v7.3');

 case 'G'

  % In this step, I compute more accurate estimates of the Bayes factors for
  % enrichment of the MHC and the xMHC using a larger range of settings for
  % the enrichment parameter, and using a more fine-grained grid for the
  % hyperparameters, with the range of the hyperparameters selected based on
  % the results of Stage F. I also output the posterior statistics about
  % individual SNPs conditioned on enrichment the MHC and xMHC.
  theta0 = (-6:0.1:-4.5)';
  theta  = (3:0.1:5)';
  
  % SPECIFY PATHWAY HYPOTHESES.
  paths = [ 3322 3323 3145 3146 3195 ];
  H     = genhmatrix(num2cell(paths),length(pathway.label));

  % Load the posterior statistics under the null hypothesis (Stage D).
  null = getnullstats('/tmp/pcarbo/pathway-t1d-D.mat',theta0);

  % COMPUTE BAYES FACTORS FOR ENRICHMENT OF GENE SETS.
  % Note that there is no need to adopt the modified analysis here (as
  % implemented in function VARPATHBIN2) because I am already analyzing SNPs
  % in the pathway---that is, SNPs in the MHC region---only after first
  % estimating posterior statistics of the SNPs outside the MHC, which were
  % obtained under the null hypothesis.
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
  [BF logw1 alpha mu s] = varpathbin(X,H,A,y,sd^2,theta0,theta,...
                                     null.logw,null.alpha,null.mu,...
                                     null.eta,gene,pathway);

  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('/tmp/pcarbo/pathway-t1d-G.mat','seed','sd','da','A','H',...
       'theta0','theta','BF','logw1','alpha','mu','s','-v7.3');

 case 'H'

  % In Stage H, I quantify support for a single enrichment hypothesis: genes
  % in the "extended" MHC that are outside the "classical" MHC region are
  % enriched at a different rate compared to genes lying within the
  % classical MHC. Since I'm only computing the Bayes factor for a single
  % enrichment hypothesis, it won't take too long to calculate an accurate
  % estimate of the Bayes factor using a fine-grained grid over the
  % hyperparameters.
  theta0   = (-6:0.1:-4.5)';
  theta    = (0:0.1:5)';
  thetaMHC = 4.6;

  % Get the posterior statistics for all SNPs given that the MHC region is
  % enriched for disease risk factors. These posterior statistics are
  % obtained by combining the posterior statistics computed in Stages D and
  % G of the analysis.
  null = getnullstats('/tmp/pcarbo/pathway-t1d-D.mat',theta0);
  mhcnull = load('/tmp/pcarbo/pathway-t1d-G.mat');
  [logw0 alpha0 mu0 s0] = getpathwaynull(null,mhcnull,theta0,thetaMHC,1);

  % Remove SNPs in the MHC region from THE pathway annotations.
  snps = find(A(:,3322));
  A(snps,:) = 0;

  % SPECIFY PATHWAY HYPOTHESES.
  m = length(pathway.label);
  H = sparse(m,1);
  H(3323,1) = 1;

  % COMPUTE BAYES FACTORS FOR ENRICHMENT OF GENE SETS.
  % Note that there is no need to adopt the modified analysis here (as
  % implemented in function VARPATHBIN2) because I am already analyzing SNPs
  % in the pathway---that is, SNPs in the MHC region---only after first
  % estimating posterior statistics of the SNPs outside the MHC, which were
  % obtained under the null hypothesis.
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
  [BF logw1] = varpathbin(X,H,A,y,sd^2,theta0,theta,logw0,alpha0,mu0,...
                          null.eta,gene,pathway);

  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('/tmp/pcarbo/pathway-t1d-H.mat','seed','sd','da','A','H',...
       'theta0','theta','thetaMHC','BF','logw1','-v7.3');

 case 'I'

  % The goal of this step of the analysis is to narrow the search for enriched
  % gene sets by obtaining rough estimates of the Bayes factors for all
  % candidate gene sets conditioned on enrichment of SNPs within the
  % MHC. Here, I compute the Bayes factors using a coarse grid for the
  % hyperparameters, and over a restricted range of the most plausible
  % hyperparameter settings. Later, I will use a more fine-grained grid over
  % a wider range of hyperparameters settings to refine my numerical
  % estimates of the Bayes factors for the most promising candidate
  % pathways.
  theta0   = (-6:0.5:-4.5)';
  theta    = (0:0.5:5)';
  thetaMHC = 4.6;

  % SPECIFY PATHWAY HYPOTHESES.
  % Remove the MHC and xMHC as candidate pathways.
  H = getinitpaths(gene,pathway);
  I = H(3322,:) | H(3323,:);
  H = H(:,find(~I));

  % Get the posterior statistics for all SNPs given that the MHC region is
  % enriched for disease risk factors. These posterior statistics are
  % obtained by combining the posterior statistics computed in Stages D and
  % G of the analysis.
  null = getnullstats('/tmp/pcarbo/pathway-t1d-D.mat',theta0);
  mhcnull = load('/tmp/pcarbo/pathway-t1d-G.mat');
  [logw0 alpha0 mu0 s0] = getpathwaynull(null,mhcnull,theta0,thetaMHC,1);

  % Load the genotype data from all SNPs except SNPs with the "extended"
  % MHC region.
  clear X
  load('/tmp/pcarbo/t1dnomhc.mat');

  % COMPUTE BAYES FACTORS FOR ENRICHMENT OF GENE SETS.
  % I remove SNPs in the "extended" MHC region from the analysis. This can be
  % seen as an approximation to the analysis in which (1) we ignore SNPs
  % within the MHC when examining SNPs outside the MHC, and (2) we have
  % separate rates of enrichment for SNPs within the MHC, and for SNPs in
  % the candidate pathway that are outside the MHC.
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
  snps       = othersnps;
  [BF logw1] = varpathbin(X,H,A(snps,:),y,sd^2,theta0,theta,logw0,...
                          alpha0(snps,:),mu0(snps,:),null.eta,...
                          gene,pathway);

  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('/tmp/pcarbo/pathway-t1d-I.mat','seed','sd','da','A','H',...
       'mhcsnps','othersnps','theta0','theta','thetaMHC','BF',...
       'logw1','-v7.3');

 case 'J'

  % In this step, I compute more accurate estimates of the Bayes factors for
  % selected gene sets using a fine-grained grid for the hyperparameters.
  theta0   = (-6:0.1:-5)';
  theta    = (2:0.1:4.5)';
  thetaMHC = 4.6;

  % SPECIFY PATHWAY HYPOTHESES.
  % I take all pathways from Stage I with Bayes factors greater than 2e8. 
  paths = [ 1250 1310 1312 1343 1354 1360 1382 1385 1390 1397 ...
            1407 1437 1757 2088 2597 2608 2653 ];
  H     = genhmatrix(num2cell(paths),length(pathway.label));

  % Get the posterior statistics for all SNPs given that the MHC region is
  % enriched for disease risk factors. These posterior statistics are
  % obtained by combining the posterior statistics computed in Stages D and
  % G of the analysis.
  null = getnullstats('/tmp/pcarbo/pathway-t1d-D.mat',theta0);
  mhcnull = load('/tmp/pcarbo/pathway-t1d-G.mat');
  [logw0 alpha0 mu0 s0] = getpathwaynull(null,mhcnull,theta0,thetaMHC,1);

  % Load the genotype data from all SNPs except SNPs with the "extended"
  % MHC region.
  clear X
  load('/tmp/pcarbo/t1dnomhc.mat');

  % COMPUTE BAYES FACTORS FOR ENRICHMENT OF GENE SETS.
  % I remove SNPs in the "extended" MHC region from the analysis. This can be
  % seen as an approximation to the analysis in which (1) we ignore SNPs
  % within the MHC when examining SNPs outside the MHC, and (2) we have
  % separate rates of enrichment for SNPs within the MHC, and for SNPs in
  % the candidate pathway that are outside the MHC.
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
  snps       = othersnps;
  [BF logw1] = varpathbin(X,H,A(snps,:),y,sd^2,theta0,theta,logw0,...
                          alpha0(snps,:),mu0(snps,:),null.eta,...
                          gene,pathway);

  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('/tmp/pcarbo/pathway-t1d-J.mat','seed','sd','da','A','H',...
       'mhcsnps','othersnps','theta0','theta','thetaMHC','BF',...
       'logw1','-v7.3');

 case 'K'

  % Next, I search for pairs of enriched pathways conditioned on enrichment of
  % MHC SNPs by computing rough estimates of the Bayes factors for candidate
  % pairs that are selected based on the results of Stage J.
  theta0   = (-6:0.2:-5.4)';
  theta    = (0:0.5:5)';
  thetaMHC = 4.6;

  % SPECIFY PATHWAY HYPOTHESES.
  % I compute Bayes factors for the top pathway from Stage J, "IL-2 mediated
  % signaling events", combined with the remaining pathways with Bayes
  % factors greater than 10 according to ranking from Stage I.
  paths = [ 12 16 18 20 21 34 55 60 63 77 91 92 93 95 103 125 126 201 222 ...
            258 263 279 287 288 295 296 299 305 306 307 312 317 318 319 ...
            320 321 322 324 327 328 335 344 354 363 365 368 369 370 371 ...
            381 384 385 386 387 402 408 440 447 455 532 533 549 550 551 ...
            556 557 562 597 599 600 602 603 604 605 606 607 608 609 612 ...
            619 621 624 627 629 630 714 725 726 748 895 913 960 963 964 ...
            965 968 969 1136 1176 1245 1250 1252 1253 1254 1265 1268 1273 ...
            1277 1281 1285 1286 1297 1310 1312 1315 1317 1322 1327 1328 ...
            1331 1340 1342 1343 1354 1355 1359 1360 1363 1365 1369 1371 ...
            1376 1381 1382 1384 1385 1390 1391 1396 1397 1398 1400 1405 ...
            1407 1412 1416 1417 1418 1419 1421 1422 1426 1427 1429 1432 ...
            1436 1437 1441 1443 1447 1448 1452 1453 1477 1490 1507 1508 ...
            1523 1525 1526 1527 1535 1566 1582 1583 1590 1592 1595 1630 ...
            1641 1644 1666 1673 1675 1688 1698 1700 1727 1728 1729 1756 ...
            1757 1758 1762 1771 1773 1777 1781 1782 1789 1791 1792 1795 ...
            1796 1805 1830 1833 1847 1887 1890 1892 1896 1910 1995 2049 ...
            2053 2055 2056 2069 2080 2086 2088 2091 2093 2106 2111 2117 ...
            2127 2128 2131 2132 2134 2138 2141 2147 2149 2154 2157 2165 ...
            2166 2184 2197 2202 2205 2208 2210 2213 2225 2226 2231 2233 ...
            2234 2237 2242 2249 2274 2291 2292 2295 2303 2304 2316 2322 ...
            2325 2326 2336 2337 2371 2388 2389 2413 2417 2418 2419 2420 ...
            2423 2431 2441 2523 2530 2559 2574 2576 2578 2584 2590 2591 ...
            2596 2597 2602 2618 2624 2636 2638 2639 2641 2642 2643 2653 ...
            2659 2665 2668 2669 2686 2688 2705 2706 2707 2708 2748 2801 ...
            2803 2827 2833 2838 2839 2840 2841 2860 2863 2873 2891 2936 ...
            2937 2938 2939 2956 2957 2960 2962 2982 2983 2985 2986 2987 ...
            2988 2989 2990 2991 2992 2993 2994 2995 3004 3023 3024 3025 ...
            3027 3028 3029 3055 3063 3065 3069 3070 3072 3073 3074 3075 ...
            3076 3079 3081 3082 3083 3085 3088 3090 3093 3094 3095 3116 ...
            3134 3136 3141 3144 3145 3146 3148 3151 3153 3158 3161 3164 ...
            3169 3174 3180 3187 3188 3189 3190 3191 3192 3198 3202 3204 ...
            3206 3214 3215 3217 3218 3219 3220 3221 3233 3240 3266 3269 3319 ];
  H         = genhmatrix(num2cell(paths),length(pathway.label));
  H(2608,:) = 1;

  % Get the posterior statistics for all SNPs given that the MHC region is
  % enriched for disease risk factors. These posterior statistics are
  % obtained by combining the posterior statistics computed in Stages D and
  % G of the analysis.
  null = getnullstats('/tmp/pcarbo/pathway-t1d-D.mat',theta0);
  mhcnull = load('/tmp/pcarbo/pathway-t1d-G.mat');
  [logw0 alpha0 mu0 s0] = getpathwaynull(null,mhcnull,theta0,thetaMHC,1);

  % Load the genotype data from all SNPs except SNPs with the "extended"
  % MHC region.
  clear X
  load('/tmp/pcarbo/t1dnomhc.mat');

  % COMPUTE BAYES FACTORS FOR ENRICHMENT OF GENE SETS.
  % I remove SNPs in the "extended" MHC region from the analysis. This can be
  % seen as an approximation to the analysis in which (1) we ignore SNPs
  % within the MHC when examining SNPs outside the MHC, and (2) we have
  % separate rates of enrichment for SNPs within the MHC, and for SNPs in
  % the candidate pathway that are outside the MHC.
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
  snps       = othersnps;
  [BF logw1] = varpathbin(X,H,A(snps,:),y,sd^2,theta0,theta,logw0,...
                          alpha0(snps,:),mu0(snps,:),null.eta,...
                          gene,pathway);

  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('/tmp/pcarbo/pathway-t1d-K.mat','seed','sd','da','A','H',...
       'mhcsnps','othersnps','theta0','theta','thetaMHC','BF',...
       'logw1','-v7.3');

 case 'L'

  % In this step, I compute more accurate estimates of Bayes factors for some
  % of the enrichment hypotheses in which pairs of pathways are enriched,
  % using a more fine-grained grid for the hyperparameters, with the range
  % of hyperparameters selected based on the results of Stage K. I also
  % output posterior statistics about individual SNPs conditioned on
  % enrichment of the selected gene sets.
  theta0   = (-6:0.1:-5.4)';
  theta    = (3:0.1:4.5)';
  thetaMHC = 4.6;

  % SPECIFY PATHWAY HYPOTHESES.
  % I recompute Bayes factors for the enrichment hypotheses with Bayes
  % factors greater than 1e13 according to the results of Stage K.
  paths     = [ 16 18 91 125 126 201 258 295 368 381 621 895 1136 1254 ...
                1285 1312 1365 1429 1630 1700 1757 1805 2111 2149 2208 ...
                2210 2325 2326 2441 2576 2597 2705 2706 2801 2838 2873 ...
                2937 3063 3079 3144 3206 3240 ];
  H         = genhmatrix(num2cell(paths),length(pathway.label));
  H(2608,:) = 1;

  % Get the posterior statistics for all SNPs given that the MHC region is
  % enriched for disease risk factors. These posterior statistics are
  % obtained by combining the posterior statistics computed in Stages D and
  % G of the analysis.
  null = getnullstats('/tmp/pcarbo/pathway-t1d-D.mat',theta0);
  mhcnull = load('/tmp/pcarbo/pathway-t1d-G.mat');
  [logw0 alpha0 mu0 s0] = getpathwaynull(null,mhcnull,theta0,thetaMHC,1);

  % Load the genotype data from all SNPs except SNPs with the "extended"
  % MHC region.
  clear X
  load('/tmp/pcarbo/t1dnomhc.mat');

  % COMPUTE BAYES FACTORS FOR ENRICHMENT OF GENE SETS.
  % I remove SNPs in the "extended" MHC region from the analysis. This can be
  % seen as an approximation to the analysis in which (1) we ignore SNPs
  % within the MHC when examining SNPs outside the MHC, and (2) we have
  % separate rates of enrichment for SNPs within the MHC, and for SNPs in
  % the candidate pathway that are outside the MHC.
  [m n] = size(H);
  fprintf('STAGE L: integrated analysis of SNPs and pathways, ');
  fprintf('given additive effects\n');
  fprintf('distribution with mean 0 and standard deviation %0.2f, ',sd);
  fprintf('and with SNP-gene\n');
  fprintf('annotation distance of %d kb. ',da/1e3);
  fprintf('Computing Bayes factors for %d enrichment\n',n);
  fprintf('hypotheses, with %d settings for theta0 ranging ',numel(theta0));
  fprintf('from %0.2f to %0.2f, and\n',min(theta0),max(theta0));
  fprintf('with %d settings for theta ranging from %0.2f to %0.2f.\n',...
          numel(theta),min(theta),max(theta));
  snps = othersnps;
  [BF logw1 alpha mu s] = ...
      varpathbin(X,H,A(snps,:),y,sd^2,theta0,theta,logw0,...
                 alpha0(snps,:),mu0(snps,:),null.eta,...
                 gene,pathway);

  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('/tmp/pcarbo/pathway-t1d-L.mat','seed','sd','da','A','H',...
       'mhcsnps','othersnps','theta0','theta','thetaMHC','BF',...
       'logw1','alpha','mu','s','-v7.3');

 case 'M'

  % In the previous step of the analysis, we found the greatest support for
  % the hypothesis that the "IL2-mediated signaling events" and "ErbB
  % receptor signaling network" are enriched. In this step, I search for
  % hypotheses in which three (3) pathways are enriched conditioned on
  % enrichment of SNPs within the MHC by computing rough estimates of the
  % Bayes factors for enrichment hypotheses in which "IL2-mediated signaling
  % events", "ErbB receptor signaling network" and a third pathway are
  % enriched.
  theta0   = (-6:0.2:-5.4)';
  theta    = (0:0.5:5)';
  thetaMHC = 4.6;

  % SPECIFY PATHWAY HYPOTHESES.
  % I compute Bayes factors for the top pair of pathways from Stage L,
  % "IL2-mediated signaling events" and "ErbB receptor signaling network",
  % combined with the remaining pathways with Bayes factors greater than 10
  % according to ranking we obtained from Stage I.
  paths = [ 12 16 18 20 21 34 55 60 63 77 91 92 93 95 103 125 126 201 222 ...
            258 263 279 287 288 295 296 299 305 306 307 312 317 318 319 ...
            320 321 322 324 327 328 335 344 354 363 365 368 369 370 371 ...
            381 384 385 386 387 402 408 440 447 455 532 533 549 550 551 ...
            556 557 562 597 599 600 602 603 604 605 606 607 608 609 612 ...
            619 621 624 627 629 630 714 725 726 748 895 913 960 963 964 ...
            965 968 969 1136 1176 1245 1250 1252 1253 1254 1265 1268 1273 ...
            1277 1281 1285 1286 1297 1310 1312 1315 1317 1322 1327 1328 ...
            1331 1340 1342 1343 1354 1355 1359 1360 1363 1365 1369 1371 ...
            1376 1381 1382 1384 1385 1390 1391 1396 1397 1398 1400 1405 ...
            1407 1412 1416 1417 1418 1419 1421 1422 1426 1427 1429 1432 ...
            1436 1437 1441 1443 1447 1448 1452 1453 1477 1490 1507 1508 ...
            1523 1525 1526 1527 1535 1566 1582 1583 1590 1592 1595 1630 ...
            1641 1644 1666 1673 1675 1688 1698 1700 1727 1728 1729 1756 ...
            1757 1758 1762 1771 1773 1777 1781 1782 1789 1791 1792 1795 ...
            1796 1805 1830 1833 1847 1887 1890 1892 1896 1910 1995 2049 ...
            2053 2055 2056 2069 2080 2086 2088 2091 2093 2106 2111 2117 ...
            2127 2128 2131 2132 2134 2138 2141 2147 2149 2154 2157 2165 ...
            2166 2184 2197 2202 2205 2208 2210 2213 2225 2226 2231 2233 ...
            2234 2237 2242 2249 2274 2291 2292 2295 2303 2304 2316 2322 ...
            2325 2326 2336 2337 2371 2388 2389 2413 2417 2418 2419 2420 ...
            2423 2431 2441 2523 2530 2559 2574 2578 2584 2590 2591 2596 ...
            2597 2602 2618 2624 2636 2638 2639 2641 2642 2643 2653 2659 ...
            2665 2668 2669 2686 2688 2705 2706 2707 2708 2748 2801 2803 ...
            2827 2833 2838 2839 2840 2841 2860 2863 2873 2891 2936 2937 ...
            2938 2939 2956 2957 2960 2962 2982 2983 2985 2986 2987 2988 ...
            2989 2990 2991 2992 2993 2994 2995 3004 3023 3024 3025 3027 ...
            3028 3029 3055 3063 3065 3069 3070 3072 3073 3074 3075 3076 ...
            3079 3081 3082 3083 3085 3088 3090 3093 3094 3095 3116 3134 ...
            3136 3141 3144 3145 3146 3148 3151 3153 3158 3161 3164 3169 ...
            3174 3180 3187 3188 3189 3190 3191 3192 3198 3202 3204 3206 ...
            3214 3215 3217 3218 3219 3220 3221 3233 3240 3266 3269 3319 ];
  H         = genhmatrix(num2cell(paths),length(pathway.label));
  H(2576,:) = 1;
  H(2608,:) = 1;

  % Get the posterior statistics for all SNPs given that the MHC region is
  % enriched for disease risk factors. These posterior statistics are
  % obtained by combining the posterior statistics computed in Stages D and
  % G of the analysis.
  null = getnullstats('/tmp/pcarbo/pathway-t1d-D.mat',theta0);
  mhcnull = load('/tmp/pcarbo/pathway-t1d-G.mat');
  [logw0 alpha0 mu0 s0] = getpathwaynull(null,mhcnull,theta0,thetaMHC,1);

  % Load the genotype data from all SNPs except SNPs with the "extended"
  % MHC region.
  clear X
  load('/tmp/pcarbo/t1dnomhc.mat');

  % COMPUTE BAYES FACTORS FOR ENRICHMENT OF GENE SETS.
  % I remove SNPs in the "extended" MHC region from the analysis. This can be
  % seen as an approximation to the analysis in which (1) we ignore SNPs
  % within the MHC when examining SNPs outside the MHC, and (2) we have
  % separate rates of enrichment for SNPs within the MHC, and for SNPs in
  % the candidate pathway that are outside the MHC.
  [m n] = size(H);
  fprintf('STAGE M: integrated analysis of SNPs and pathways, ');
  fprintf('given additive effects\n');
  fprintf('distribution with mean 0 and standard deviation %0.2f, ',sd);
  fprintf('and with SNP-gene\n');
  fprintf('annotation distance of %d kb. ',da/1e3);
  fprintf('Computing Bayes factors for %d enrichment\n',n);
  fprintf('hypotheses, with %d settings for theta0 ranging ',numel(theta0));
  fprintf('from %0.2f to %0.2f, and\n',min(theta0),max(theta0));
  fprintf('with %d settings for theta ranging from %0.2f to %0.2f.\n',...
          numel(theta),min(theta),max(theta));
  snps       = othersnps;
  [BF logw1] = varpathbin(X,H,A(snps,:),y,sd^2,theta0,theta,logw0,...
                          alpha0(snps,:),mu0(snps,:),null.eta,...
                          gene,pathway);

  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('/tmp/pcarbo/pathway-t1d-M.mat','seed','sd','da','A','H',...
       'mhcsnps','othersnps','theta0','theta','thetaMHC','BF',...
       'logw1','-v7.3');

 case 'N'

  % In this step, I compute more accurate estimates of Bayes factors for
  % some of the enrichment hypotheses in which three (3) pathways are
  % enriched using a more fine-grained grid for the hyperparameters, with
  % the range of hyperparameters selected based on the outcome of Stage M.
  theta0   = (-6:0.1:-5.4)';
  theta    = (3:0.1:4.5)';
  thetaMHC = 4.6;

  % SPECIFY PATHWAY HYPOTHESES.
  % I recompute Bayes factors for the 8 enrichment hypotheses with Bayes
  % factors greater than 1e16 according to the results of Stage M.
  paths     = [ 125 381 1136 1254 1312 1700 2326 2597 ];
  H         = genhmatrix(num2cell(paths),length(pathway.label));
  H(2576,:) = 1;
  H(2608,:) = 1;

  % Get the posterior statistics for all SNPs given that the MHC region is
  % enriched for disease risk factors. These posterior statistics are
  % obtained by combining the posterior statistics computed in Stages D and
  % G of the analysis.
  null = getnullstats('/tmp/pcarbo/pathway-t1d-D.mat',theta0);
  mhcnull = load('/tmp/pcarbo/pathway-t1d-G.mat');
  [logw0 alpha0 mu0 s0] = getpathwaynull(null,mhcnull,theta0,thetaMHC,1);

  % Load the genotype data from all SNPs except SNPs with the "extended"
  % MHC region.
  clear X
  load('/tmp/pcarbo/t1dnomhc.mat');

  % COMPUTE BAYES FACTORS FOR ENRICHMENT OF GENE SETS.
  % I remove SNPs in the "extended" MHC region from the analysis. This can be
  % seen as an approximation to the analysis in which (1) we ignore SNPs
  % within the MHC when examining SNPs outside the MHC, and (2) we have
  % separate rates of enrichment for SNPs within the MHC, and for SNPs in
  % the candidate pathway that are outside the MHC.
  [m n] = size(H);
  fprintf('STAGE N: integrated analysis of SNPs and pathways, ');
  fprintf('given additive effects\n');
  fprintf('distribution with mean 0 and standard deviation %0.2f, ',sd);
  fprintf('and with SNP-gene\n');
  fprintf('annotation distance of %d kb. ',da/1e3);
  fprintf('Computing Bayes factors for %d enrichment\n',n);
  fprintf('hypotheses, with %d settings for theta0 ranging ',numel(theta0));
  fprintf('from %0.2f to %0.2f, and\n',min(theta0),max(theta0));
  fprintf('with %d settings for theta ranging from %0.2f to %0.2f.\n',...
          numel(theta),min(theta),max(theta));
  snps       = othersnps;
  [BF logw1] = varpathbin(X,H,A(snps,:),y,sd^2,theta0,theta,logw0,...
                          alpha0(snps,:),mu0(snps,:),null.eta,...
                          gene,pathway);

  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('/tmp/pcarbo/pathway-t1d-N.mat','seed','sd','da','A','H',...
       'mhcsnps','othersnps','theta0','theta','thetaMHC','BF',...
       'logw1','-v7.3');

 case 'O'

  % In this step, I compute posterior statistics about individual SNPs
  % conditioned on selected enrichment hypotheses in which one pathway is
  % enriched.
  theta0   = (-5.8:0.1:-5.4)';
  theta    = (3:0.1:4.5)';
  thetaMHC = 4.6;

  % SPECIFY PATHWAY HYPOTHESES.
  % I take all pathways from Stage J with Bayes factors greater than 1e10.
  paths = [ 1312 1343 1354 1360 1390 2088 2597 2608 ];
  H     = genhmatrix(num2cell(paths),length(pathway.label));

  % Get the posterior statistics for all SNPs given that the MHC region is
  % enriched for disease risk factors. These posterior statistics are
  % obtained by combining the posterior statistics computed in Stages D and
  % G of the analysis.
  null = getnullstats('/tmp/pcarbo/pathway-t1d-D.mat',theta0);
  mhcnull = load('/tmp/pcarbo/pathway-t1d-G.mat');
  [logw0 alpha0 mu0 s0] = getpathwaynull(null,mhcnull,theta0,thetaMHC,1);

  % Load the genotype data from all SNPs except SNPs with the "extended"
  % MHC region.
  clear X
  load('/tmp/pcarbo/t1dnomhc.mat');

  % COMPUTE BAYES FACTORS FOR ENRICHMENT OF GENE SETS.
  % I remove SNPs in the "extended" MHC region from the analysis. This can be
  % seen as an approximation to the analysis in which (1) we ignore SNPs
  % within the MHC when examining SNPs outside the MHC, and (2) we have
  % separate rates of enrichment for SNPs within the MHC, and for SNPs in
  % the candidate pathway that are outside the MHC.
  [m n] = size(H);
  fprintf('STAGE O: integrated analysis of SNPs and pathways, ');
  fprintf('given additive effects\n');
  fprintf('distribution with mean 0 and standard deviation %0.2f, ',sd);
  fprintf('and with SNP-gene\n');
  fprintf('annotation distance of %d kb. ',da/1e3);
  fprintf('Computing Bayes factors for %d enrichment\n',n);
  fprintf('hypotheses, with %d settings for theta0 ranging ',numel(theta0));
  fprintf('from %0.2f to %0.2f, and\n',min(theta0),max(theta0));
  fprintf('with %d settings for theta ranging from %0.2f to %0.2f.\n',...
          numel(theta),min(theta),max(theta));
  snps = othersnps;
  [BF logw1 alpha mu s] = ...
      varpathbin(X,H,A(snps,:),y,sd^2,theta0,theta,logw0,...
                 alpha0(snps,:),mu0(snps,:),null.eta,...
                 gene,pathway);

  % SAVE RESULTS.
  fprintf('Saving results.\n');
  save('/tmp/pcarbo/pathway-t1d-O.mat','seed','sd','da','A','H',...
       'mhcsnps','othersnps','theta0','theta','thetaMHC','BF',...
       'logw1','alpha','mu','s','-v7.3');
end

