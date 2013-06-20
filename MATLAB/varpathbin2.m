% [BF,LOGW1,ALPHA,MU,S] = VARPATHBIN2(X,H,A,Y,S,SA,THETA0,THETA,LOGW0,...
% ALPHA0,MU0,ETA0,GENE,PATHWAY) does the same thing as VARPATHBIN, except
% that we compute posterior statistics for SNPs annotated by S (i.e. SNPs i
% for which S(i) == 1) ignoring SNPs that are not annotated by S.
function [BF, logw1, alpha, mu, s] = ...
        varpathbin2 (X, H, A, y, S, sa, theta0, theta, logw0, ...
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
    a     = spones(sumcols(A(:,paths)));
    snps  = find(a);
    snps1 = find(a & S);
    snps2 = find(a & ~S);
    genes = find(sumcols(pathway.genes(:,paths)) & inref(gene));

    % Display some information about the enrichment hypothesis.
    fprintf('ENRICHMENT HYPOTHESIS %d out of %d ',i,n);
    fprintf('implicates %d genes and\n',length(genes));
    fprintf('%d + %d = %d SNPs:\n',length(snps1),...
            length(snps2),length(snps));
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
    if length(snps) == 0
      BF(i) = 0;
    else

      % Compute the Bayes factor.
      if length(snps1) == 0 | length(snps2) == 0
        [BF(i) logw1(:,:,i) alphapath mupath spath] = ...
            bayesfactorbin(X(:,snps),y,sa,theta0,theta,...
                           logw0,alpha0(snps,:),mu0(snps,:),eta0,Xr0);
      else
        [BF(i) logw1(:,:,i) alphapath mupath spath] = ...
            bayesfactorbin2(X(:,snps),y,find(S(snps)),find(~S(snps)),...
                            sa,theta0,theta,logw0,alpha0(snps,:),...
                            mu0(snps,:),eta0,Xr0);
      end

      % Store posterior statistics about individual SNPs, if requested.
      if nargout > 2
        alpha{i}(snps,:) = reshape(alphapath,length(snps),n0*n1);
        mu{i}(snps,:)    = reshape(mupath,length(snps),n0*n1);
        s{i}(snps,:)     = reshape(spath,length(snps),n0*n1);
      end
    end
    fprintf('BF = %0.2e\n',BF(i));      
  end
