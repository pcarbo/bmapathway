% [BF,LOGW1,ALPHA,MU,S] = BAYESFACTORBIN2(X,Y,A,B,SA,THETA0,THETA,LOGW0,...
% ALPHA0,MU0,ETA0,XR0) does the same thing as BAYESFACTORBIN, except that
% posterior statistics for the first set of SNPs (A) are estimated by
% ignoring the second set of SNPs (B).
function [BF, logw1, alpha, mu, s] = ...
      bayesfactorbin2 (X, y, A, B, sa, theta0, theta, ...
                       logw0, alpha0, mu0, eta0, Xr0)

  % Get the number of samples (n), the number of SNPs assigned to the enriched
  % pathway (p), the number of settings of the genome-wide log-odds (n0),
  % and the number of settings of the enrichment parameter (n1).
  [n p] = size(X);
  n0    = numel(theta0);
  n1    = numel(theta);

  % Set a random initialization of the variational parameters for each
  % combination of the hyperparameters.
  alpha = rand(p,n0,n1);
  alpha = alpha ./ repmat(sum(alpha),p,1);
  mu    = randn(p,n0,n1);

  % First get the best initialization for the variational parameters.
  fprintf('Finding best initialization for %d combinations ',n0*n1);
  fprintf('of hyperparameters.\n');
  [logw1 alpha mu s] = ...
      outerloopbayesfactorbin2(X,y,A,B,sa,theta0,theta,alpha,mu,...
                               logw0,alpha0,mu0,eta0,Xr0);

  % Choose an initialization common to all the runs of the coordinate
  % ascent algorithm. This is chosen from the hyperparameters with the
  % highest variational estimate of the importance weight.
  [ans i] = max(logw1(:));
  alpha   = repmat(alpha(:,i),[1 n0 n1]);
  mu      = repmat(mu(:,i),[1 n0 n1]);

  % Compute the unnormalized log-importance weights.
  fprintf('Computing importance weights for %d combinations ',n0*n1);
  fprintf('of hyperparameters.\n');
  [logw1 alpha mu s] = ...
      outerloopbayesfactorbin2(X,y,A,B,sa,theta0,theta,alpha,mu,...
                               logw0,alpha0,mu0,eta0,Xr0);

  % COMPUTE BAYES FACTOR. 
  % Compute the marginal log-likelihood under the null hypothesis using
  % importance sampling.
  c     = max(logw0(:));
  logZ0 = c + log(mean(exp(logw0(:) - c)));

  % Compute the marginal log-likelihood under the enrichment hypothesis using
  % importance sampling.
  c     = max(logw1(:));
  logZ1 = c + log(mean(exp(logw1(:) - c)));

  % Get the numerical estimate of the Bayes factor.
  BF = exp(logZ1 - logZ0);

