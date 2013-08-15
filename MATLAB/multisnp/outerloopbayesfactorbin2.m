% [LOGW1,ALPHA,MU,S] = OUTERLOOPBAYESFACTORBIN2(X,Y,A,B,SA,THETA0,THETA,...
% ALPHA,MU,LOGW0,ALPHA0,MU0,ETA0,XR0) does the same thing as
% OUTERLOOPBAYESFACTORBIN, except that posterior statistics for the first
% set of SNPs (A) are estimated by ignoring the second set of SNPs (B).
function [logw1, alpha, mu, s] = ...
    outerloopbayesfactorbin2 (X, y, A, B, sa, theta0, theta, ...
                              alpha, mu, logw0, alpha0, mu0, eta0, Xr0)

  % Get the number of samples (n), the number of SNPs assigned to the enriched
  % pathway (p), the number of settings of the genome-wide log-odds (n0),
  % and the number of settings of the enrichment parameter (n1).
  [n p] = size(X);
  n0    = numel(theta0);
  n1    = numel(theta);

  % Initialize storage for log-importance weights under the alternative
  % hypothesis, and for the variances of the additive effects.
  logw1 = zeros(n0,n1);
  s     = zeros(p,n0,n1);

  % Repeat for each setting of the genome-wide log-odds (theta0).
  iter = 0;
  for i = 1:n0

    % Compute YHAT1 = Y - UHAT*X(:,S)*R(S), the phenotype samples
    % "adjusted" for additive effects of SNPs that are outside the
    % pathway. Here I've defined S to be the set of SNPs outside the
    % pathway, and UHAT is defined in the Bayesian Analysis (2012) paper;
    % it is a function of the variational parameters ETA.
    u     = slope(eta0(:,i));
    Xrp   = double(X*(alpha0(:,i).*mu0(:,i)));
    Xrs   = Xr0(:,i) - Xrp;    
    yhat1 = y - (u.*Xrs - u*(u'*Xrs)/sum(u));

    % Compute the variational lower bound to the marginal log-likelihood under
    % the null hypothesis when there is no enrichment, only for SNPs in the
    % pathway. Note that VARBVSBIN defines the log-odds ratio using the
    % natural logarithm, so we need to multiply THETA0 by LOG(10).
    logodds = log(10) * theta0(i);
    options = struct('alpha',alpha0(:,i),'mu',mu0(:,i),'eta',eta0(:,i),...
		     'fixed_eta',true,'verbose',false);
    F0 = varbvsbin(X,yhat1,sa,logodds,options);  

    % Repeat for each setting of the enrichment parameter (theta).
    for j = 1:n1

      % Display the current status of the inference procedure.
      iter = iter + 1;
      fprintf('(%04d) theta0 = %+0.2f, theta = %0.2f',...
              iter,theta0(i),theta(j));
      fprintf(repmat('\b',1,35));

      % Compute the marginal log-likelihood under the alternative
      % hypothesis only for SNPs in the pathway. Note that VARBVSBIN
      % defines the log-odds ratio using the natural logarithm, so we
      % need to multiply by LOG(10).
      % 
      % Start by estimating posterior probabilities for the first set of 
      % SNPs (set A).
      logodds = log(10) * (theta0(i) + theta(j));
      options = struct('alpha',alpha(A,i,j),'mu',mu(A,i,j),'eta',...
                       eta0(:,i),'fixed_eta',true,'verbose',false);
      [lnZ alpha(A,i,j) mu(A,i,j) s(A,i,j)] = ...
	  varbvsbin(X(:,A),yhat1,sa,logodds,options);

      % Compute YHAT2 = Y - UHAT*X(:,W)*R(W), the phenotype samples
      % "adjusted" for additive effects of SNPs that are in set W =
      % UNION(A,S), where A is the first set of SNPs analyzed, and S is
      % the set of SNPs outside the pathway. UHAT is defined in the
      % Bayesian Analysis (2012) paper; it is a function of the
      % variational parameters ETA. 
      Xra   = double(X(:,A) * (alpha(A,i).*mu(A,i)));
      Xrw   = Xra + Xrs;
      yhat2 = y - (u.*Xrw - u*(u'*Xrw)/sum(u));

      % Estimate the posterior probabilities for the second set of SNPs
      % (set B) conditioned on the additive effects of the first set of
      % SNPs (set A), in addition to the SNPs outside the pathway.
      logodds = log(10) * (theta0(i) + theta(j));
      options = struct('alpha',alpha(B,i,j),'mu',mu(B,i,j),'eta',...
                       eta0(:,i),'fixed_eta',true,'verbose',false);
      [lnZ alpha(B,i,j) mu(B,i,j) s(B,i,j)] = ...
	  varbvsbin(X(:,B),yhat2,sa,logodds,options);

      % Get the variational estimate of the marginal log-likelihood under
      % the alternative hypothesis only for SNPs in the pathway.
      Xr    = double(X*(alpha(:,i,j).*mu(:,i,j)));
      stats = updatestats(X,yhat1,eta0(:,i));
      F1    = intlogit(yhat1,stats,alpha(:,i,j),mu(:,i,j),s(:,i,j),...
                       Xr,eta0(:,i)) ...
              + intgamma(logodds,alpha(:,i,j)) ...
              + intklbeta(alpha(:,i,j),mu(:,i,j),s(:,i,j),sa);    

      % Compute the importance weight under the alternative. The prior and
      % proposal for the genome-wide log-odds and enrichment parameters are
      % assumed to be uniform, so they cancel out from the importance
      % weights. For an explanation why we can decompose the variational
      % lower bound in this way, see the supporting material in the PLoS
      % Genetics manuscript.
      logw1(i,j) = logw0(i) + F1 - F0;
    end
  end
  fprintf('\n');
