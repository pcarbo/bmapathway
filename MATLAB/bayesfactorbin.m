% [BF,LOGW1] = BAYESFACTORBIN(X,Y,SA,THETA0,THETA,LOGW0,ALPHA0,...
% MU0,ETA0,XR0) computes the Bayes factor for a pathway hypothesis, assuming
% a binary trait.
%
% Input X is the genotype data for SNPs assigned to the pathway. It is an N
% x P matrix, where N is the number of samples (individuals) and P is the
% number of genetic loci, or SNPs, assigned to the enriched pathway, or
% pathways. Y is the vector of observations about the binary trait. X and Y
% should not be centered. Instead, we will account for the intercept as we
% update the variational approximation.
%
% SA is a scalar specifying the prior variance of the additive effects,
% THETA0 is an array of settings of the genome-wide log-odds, and THETA is
% an array of settings of the enrichment. For the alternative hypothesis
% that a pathway is enriched for the trait, importance weights are computed
% for all combinations of THETA0 and THETA. These importance weights are
% calculated under the assumption that the prior and proposal distributions
% for THETA0 and THETA are uniform, so they cancel out in the expression for
% the importance weights.
%
% Inputs ALPHA0, MU0, W0 and XR0 are previously calculated quantities under
% the null hypothesis that no pathways are enriched for the trait. ALPHA0
% and MU0 are variational estimates of the posterior inclusion probabilities
% and posterior means of the additive effects for each setting. They are
% each arrays of dimension P x NUMEL(THETA0). LOGW0 is an array of the same
% size as THETA0 containing the value of the log-importance weight for
% settings of THETA0 under the null hypothesis. XR0 is an N x NUMEL(THETA0)
% matrix, in which XR0(:,i) = X*(ALPHA0(:,i).*MU0(:,i)), where X is the
% genotype matrix for *all* SNPs (not just SNPs assigned to the pathway),
% and matrices ALPHA0 and MU0 give variational estimates of posterior
% expectations for all available SNPs genome-wide. ETA0 is the set of
% parameters specifying the variational approximation to the nonlinear
% logistic regression factors for each setting of THETA0. It is an N x
% NUMEL(THETA0) matrix.
%
% Output BF is the Bayes factor. To calculate this Bayes factor, we
% integrate over the hyperparameters using importance sampling under the
% assumption that SNPs outside the enriched pathways are unaffected by
% pathways (a posteriori).
%
% Output LOGW1 is the array of log-importance weights for settings of the
% hyperparameters, under the pathway hypothesis. It is a matrix of dimension
% NUMEL(THETA0) x NUMEL(THETA). 
% 
% [BF,LOGW1,ALPHA,MU,S] = BAYESFACTORBIN(...) returns posterior statistics
% about the individual SNPs in the pathway; ALPHA, MU and S are variational
% estimates of the posterior inclusion probabilities, posterior means and
% posterior variances of the additive effects given the pathway hypothesis,
% respectively for each setting of the hyperparameters. Each of these
% outputs are arrays of dimension P x NUMEL(THETA0) x NUMEL(THETA).
function [BF, logw1, alpha, mu, s] = ...
      bayesfactorbin (X, y, sa, theta0, theta, ...
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
      outerloopbayesfactorbin(X,y,sa,theta0,theta,alpha,mu,...
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
      outerloopbayesfactorbin(X,y,sa,theta0,theta,alpha,mu,...
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

