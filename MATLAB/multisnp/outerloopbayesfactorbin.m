% [LOGW1, ALPHA, MU, S] = OUTERLOOPBAYESFACTORBIN(X,Y,SA,THETA0,THETA,...
% ALPHA,MU,LOGW0,ALPHA0,MU0,ETA0,XR0) computes the unnormalized
% log-importance weights for the hyperparameters. It is used by
% BAYESFACTORBIN to implement the "outer loop" of the variational inference
% algorithm.
function [logw1, alpha, mu, s] = ...
    outerloopbayesfactorbin (X, y, sa, theta0, theta, alpha, mu, ...
                             logw0, alpha0, mu0, eta0, Xr0)

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

    % Compute YHAT = Y - UHAT*X(:,B)*R(B), the phenotype samples "adjusted" for
    % additive effects of SNPs that are outside the pathway, where B is
    % defined to be the set of SNPs outside the pathway. UHAT is defined in
    % the Bayesian Analysis (2012) paper; it is a function of the
    % variational parameters ETA.
    u    = slope(eta0(:,i));
    Xra  = double(X*(alpha0(:,i).*mu0(:,i)));
    Xrb  = Xr0(:,i) - Xra;    
    yhat = y - (u.*Xrb - u*(u'*Xrb)/sum(u));

    % Compute the variational lower bound to the marginal log-likelihood under
    % the null hypothesis when there is no enrichment, only for SNPs in the
    % pathway. Note that VARBVSBIN defines the log-odds ratio using the
    % natural logarithm, so we need to multiply THETA0 by LOG(10).
    logodds = log(10) * theta0(i);
    options = struct('alpha',alpha0(:,i),'mu',mu0(:,i),'eta',eta0(:,i),...
		     'fixed_eta',true,'verbose',false);
    F0 = varbvsbin(X,yhat,sa,logodds,options);  

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
      options = struct('alpha',alpha(:,i,j),'mu',mu(:,i,j),'eta',...
                       eta0(:,i),'fixed_eta',true,'verbose',false);
      logodds = log(10) * (theta0(i) + theta(j));
      [F1 alpha(:,i,j) mu(:,i,j) s(:,i,j)] = ...
	  varbvsbin(X,yhat,sa,logodds,options);

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
