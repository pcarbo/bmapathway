% [BF,LOGW1] = VARPATHBIN(X,H,A,Y,SA,THETA0,THETA,LOGW0,ALPHA0,MU0,ETA0,...
% GENE,PATHWAY) computes Bayes factors for pathway enrichment hypotheses
% corresponding to a binary trait.
%
% Input X is the genotype data. It is an S x P matrix, where S is the number
% of samples (individuals), and P is the number of genetic loci, or SNPs. Y
% is the vector of observations about the binary trait. It is a vector of
% length S. X and Y should not be centered. Instead, we will account for the
% intercept as we update the variational approximation.
%
% Input A is the P x M matrix of pathway annotations, where P is the number
% of SNPs and M is the number of pathways. A(i,j) = 1 means SNP i is
% assigned to pathway j; otherwise, A(i,j) = 0. Typically A is sparse.
% Input H specifies the set of enrichment hypotheses. (Which does not
% include the null hypothesis.) H is an M x N matrix, where M is the number
% of pathways, and N is the number of enrichment hypotheses, such that
% H(j,k) = 1 if and only if pathway j is enriched for associations with the
% binary trait in enrichment hypothesis k; otherwise, H(j,k) = 0.
%
% SA is the prior variance of the additive effects, THETA0 is an array of
% settings of the genome-wide log-odds, and THETA is an array of settings of
% the enrichment parameter. For the alternative hypothesis that a pathway is
% enriched for associations with the binary trait, importance weights are
% computed for all combinations of THETA0 and THETA. These importance
% weights are calculated under the assumption that the prior and proposal
% distributions are uniform for both hyperparameters THETA0 and THETA. Note
% that a residual variance parameter is not needed to model a binary trait.
%
% Inputs LOGW0, ALPHA0, MU0 and ETA0 are quantities computed under the null
% hypothesis in which no pathways are enriched for the trait. These
% quantities should be computed using MULTISNPBINHYPER.
%
% This function call returns two outputs, BF and LOGW1. BF(k) is the Bayes
% factor for enrichment hypothesis k. LOGW1(:,:,k) gives the log-importance
% weights for the hyperparameters given enrichment hypothesis k. LOGW1 is an
% array of dimension NUMEL(THETA0) x NUMEL(THETA) x N.
%
% [BF,W1,ALPHA,MU,S] = VARPATHBIN(...) returns additional statistics about
% the individual SNPs in cell arrays ALPHA, MU and S. Entries ALPHA{i}(j,k),
% MU{i}(j,k) and S{i}(j,k) are the posterior inclusion probability,
% posterior mean and posterior variance of the additive effect for SNP j
% given pathway enrichment hypothesis i and hyperparameter setting k, given
% that SNP j is assigned to at least one of the enriched pathways in
% hypothesis k; if SNP j is not assigned to any one of the enriched
% pathways, the respective entries are set to 0.
function [BF, logw1, alpha, mu, s] = ...
        varpathbin (X, H, A, y, sa, theta0, theta, logw0, ...
                    alpha0, mu0, eta0, gene, pathway)

  % Get the number of SNPs genotyped (p), the number of pathways (m), the
  % number of pathway enrichment hypotheses (n), the number settings of the
  % genome-wide log-odds (n0), and the number of settings of the enrichment
  % parameter (n1).
  [p m] = size(A);
  [m n] = size(H);
  n0    = numel(theta0);
  n1    = numel(theta);

  % Compute the matrix-vector product X*r for each setting of the
  % genome-wide log-odds.
  Xr0 = double(X*(alpha0.*mu0));

  % Initialize storage for the Bayes factors (BF) and the log-importance
  % weights for the hyperparameters given each pathway enrichment hypothesis
  % (logw1).
  BF    = zeros(n,1);
  logw1 = zeros(n0,n1,n);

  % Initialize storage for the variational estimates of the posterior
  % statistics for each setting of the hyperparameters, and for each pathway
  % hypothesis, but only if these outputs are requested. Even though it
  % isn't necessary to initialize storage at this stage, it is useful to do
  % so to make sure that we have enough space in memory to store the results
  % from all the pathway enrichment hypotheses.
  if nargout > 2
    alpha = cell(n,1);
    mu    = cell(n,1);
    s     = cell(n,1);

    % Repeat for each pathway enrichment hypothesis.
    for i = 1:n

      % Get the number of SNPs assigned to the enriched pathway(s).
      p1 = full(sum(spones(A * H(:,i))));

      % Allocate storage for the variational estimates of the posterior
      % expectations.
      M        = spalloc(p,n0*n1,p1*n0*n1);
      alpha{i} = M;
      mu{i}    = M;
      s{i}     = M;
    end
  end

  % Compute the Bayes factor for each enrichment hypothesis.
  for i = 1:n

    % Get the set of genes and the set of SNPs assigned to the enriched
    % pathway(s).
    paths = find(H(:,i));
    snps  = find(sumcols(A(:,paths)));
    genes = find(sumcols(pathway.genes(:,paths)) & inref(gene));

    % Get the number of SNPs assigned to the enriched pathways.
    p1 = length(snps);

    % Display some information about the enrichment hypothesis.
    fprintf('ENRICHMENT HYPOTHESIS %d out of %d ',i,n);
    fprintf('implicates %d genes and %d SNPs:\n',length(genes),p1);
    for t = 1:length(paths)
      j   = paths(t);
      str = pathway.label{j};
      str = str(1:min(60,length(str)));
      if strcmpi(pathway.source{j},pathway.database{j})
	fprintf('%d. %s (%s)\n',t,str,pathway.source{j});
      else
	fprintf('%d. %s (%s/%s)\n',t,str,pathway.source{j},...
                pathway.database{j});
      end
    end

    % If no SNPs are assigned to the pathway, set the Bayes factor to 0.
    if p1 == 0
      BF(i) = 0;
    else

      % Compute the Bayes factor.
      [BF(i) logw1(:,:,i) alphapath mupath spath] = ...
          bayesfactorbin(X(:,snps),y,sa,theta0,theta,logw0,...
                         alpha0(snps,:),mu0(snps,:),eta0,Xr0);

      % Store posterior statistics about individual SNPs, if requested.
      if nargout > 2
        alpha{i}(snps,:) = reshape(alphapath,p1,n0*n1);
        mu{i}(snps,:)    = reshape(mupath,p1,n0*n1);
        s{i}(snps,:)     = reshape(spath,p1,n0*n1);
      end
    end
    fprintf('BF = %0.2e\n',BF(i));      
  end
