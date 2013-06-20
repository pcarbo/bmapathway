% [LOGW,ALPHA,MU,S] = MULTISNPHYPER(X,Y,LOG10SIGMA,H,LOG10ODDS,MU0,S0) runs
% the full variational inference procedure for Bayesian variable selection
% in linear regression (for a quantitative trait). The first part to the
% algorithm searches for a good initialization of the variational
% parameters.
% 
% This inference procedure involves an inner loop and an outer loop. The
% inner loop consists of running a coordinate ascent algorithm to tighten
% the variational lower bound given a setting of the hyperparameters. The
% outer loop computes importance weights for all combinations of the
% hyperparameters.
%
% Input X is the genotype data. It is an N x P matrix, where N is the number
% of samples (individuals), and P is the number of variables (genetic loci,
% or SNPs). Y is the vector of quantitative trait data; it is a vector of
% length N. Crucially, this algorithm will only work correctly if Y and X
% are centered so that vector Y and each column of X has a mean of zero.
%
% Inputs LOG10SIGMA, H and LOG10ODDS are the hyperparameter settings.
% LOG10SIGMA is the (base 10) logarithm of the residual variance. H is the
% prior estimate of the proportion of variance explained.  LOG10ODDS is the
% (base 10) logarithm of the prior odds for inclusion; it is equal to
% LOG10ODDS = LOG10(Q./(1-Q)), where Q is the prior probability that each
% SNP is included in the linear model of Y. These three inputs must be
% arrays of the same size. For each combination of the hyperparameters, we
% compute the unnormalized log-importance weight, and store the result in
% output LOGW.
%
% We assume a uniform prior for hyperparameters LOG10SIGMA and H. The prior
% for LOG10ODDS is uniform over the interval given by the candidate values
% of LOG10ODDS. Assuming a uniform proposal distribution for all the
% hyperparameters, the importance sampling procedure will recover the
% correct distribution of H and LOG10ODDS as the number of Monte Carlo
% samples becomes large, and assuming the variational approximation closely
% matches the exact posterior. See PVE2SA for how to obtain the prior
% variance of the additive effects.
%
% Outputs ALPHA, MU and S are variational estimates of the posterior
% inclusion probabilities, and posterior means and variances of the additive
% effects (given that the variable is included in the model) for each etting
% of the hyperparameters. ALPHA, MU and S are each a matrix of dimension P x
% NS, where NS is the number of hyperparameter settings.
function [logw, alpha, mu, s] = multisnphyper (X, y, log10sigma, h, log10odds)
  
  % Get the number of participants in the study (n), the number of SNPs
  % genotyped (p), and the number of combinations of the hyperparameters
  % (ns). 
  [n p] = size(X);
  ns    = numel(h);
  
  % Set a random initialization of the variational parameters for each
  % combination of the hyperparameters.
  alpha = rand(p,ns);
  alpha = alpha ./ repmat(sum(alpha),p,1);
  mu    = randn(p,ns);

  % First get the best initialization for the variational parameters.
  fprintf('Finding best initialization for %d combinations ',ns);
  fprintf('of hyperparameters.\n');
  [logw alpha mu] = outerloophyper(X,y,alpha,mu,log10sigma,h,log10odds);
  
  % Choose an initialization common to all the runs of the coordinate ascent
  % algorithm. This is chosen from the hyperparameters with the highest
  % variational estimate of the posterior probability.
  [ans i] = max(logw(:));
  alpha   = repmat(alpha(:,i),1,ns);
  mu      = repmat(mu(:,i),1,ns);

  % Compute the unnormalized log-importance weights.
  fprintf('Computing importance weights for %d combinations ',ns);
  fprintf('of hyperparameters.\n');
  [logw alpha mu s] = outerloophyper(X,y,alpha,mu,log10sigma,h,log10odds);
