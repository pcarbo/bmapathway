function [logw, alpha, mu, s, eta] = ...
      outerloophyperbin2 (X, y, A, B, alpha, mu, eta, h, log10odds)

  % Get the number of participants in the study (n), the number of SNPs
  % genotyped (p), and the number of combinations of the hyperparameters
  % (ns).
  [n p] = size(X);
  ns    = numel(h);

  % Get the settings for the prior variance of the additive effects.
  sx = sum(var1(X));
  sa = pve2sa(sx,h,log10odds);

  % Initialize storage for the unnormalized log-importance weights, and
  % variances of the additive effects.
  logw = zeros(size(h));
  s    = zeros(p,ns);

  % Repeat for each combination of the hyperparameters.
  for i = 1:ns
    fprintf('(%03d) h = %0.3f, log10odds = %+0.2f (sd = %0.3f)\n',...
	    i,h(i),log10odds(i),sqrt(sa(i)));
  
    % Compute the unnormalized log-importance weight given values for the
    % hyperparameters, H and LOG10ODDS. Implicitly, the importance weight
    % includes these terms: the likelihood, the prior, and the proposal. The
    % proposal and prior cancel out from the expression for the importance
    % weight because both are assumed to be uniform for all the
    % hyperparameters.
    % 
    % Start by estimating posterior probabilities for the first set of
    % SNPs (set A).
    logodds = log(10) * log10odds(i);
    options = struct('alpha',alpha(A,i),'mu',mu(A,i),'eta',eta(:,i),...
                     'fixed_eta',true,'verbose',true);
    [lnZ alpha(A,i) mu(A,i) s(A,i)] = ...
	varbvsbin(X(:,A),y,sa(i),logodds,options);

    % Compute YHAT = Y - UHAT*X(:,A)*R(A), the phenotype samples
    % "adjusted" for additive effects of SNPs that are in the first set
    % (set A). UHAT is defined in the Bayesian Analysis (2012) paper; it
    % is a function of the variational parameters ETA.
    u    = slope(eta(:,i));
    Xra  = double(X(:,A) * (alpha(A,i).*mu(A,i)));
    yhat = y - (u.*Xra - u*(u'*Xra)/sum(u));

    % Estimate the posterior probabilities for the second set of SNPs (set B)
    % conditioned on the additive effects of the first set of SNPs (set A).
    options = struct('alpha',alpha(B,i),'mu',mu(B,i),'eta',eta(:,i),...
                     'fixed_eta',true,'verbose',false);
    [lnZ alpha(B,i) mu(B,i) s(B,i)] = ...
	varbvsbin(X(:,B),yhat,sa(i),logodds,options);

    % Get the variational estimate of the marginal log-likelihood.
    Xr      = double(X*(alpha(:,i).*mu(:,i)));
    stats   = updatestats(X,y,eta(:,i));
    logw(i) = intlogit(y,stats,alpha(:,i),mu(:,i),s(:,i),Xr,eta(:,i)) ...
              + intgamma(logodds,alpha(:,i)) ...
              + intklbeta(alpha(:,i),mu(:,i),s(:,i),sa(i));    

    fprintf('\n');
  end
