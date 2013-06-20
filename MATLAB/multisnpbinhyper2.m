function [logw, alpha, mu, s] = ...
        multisnpbinhyper2 (X, y, A, B, h, log10odds, eta)

  % Get the number of participants in the study (n), the number of SNPs
  % genotyped (p), and the number of combinations of the hyperparameters
  % (ns). 
  [n p] = size(X);
  ns    = numel(h);

  % Set a random initialization of the variational parameters for each
  % combination of the hyperparameters, or use the initialization
  % provided by the input arguments.
  alpha = rand(p,ns);
  alpha = alpha ./ repmat(sum(alpha),p,1);
  mu    = randn(p,ns);

  % First get the best initialization for the variational parameters.
  fprintf('Finding best initialization for %d combinations ',ns);
  fprintf('of hyperparameters.\n');
  [logw alpha mu s] = outerloophyperbin2(X,y,A,B,alpha,mu,eta,h,log10odds);

  % Choose an initialization common to all the runs of the coordinate ascent
  % algorithm. This is chosen from the hyperparameters with the highest
  % variational estimate of the posterior probability.
  [ans i] = max(logw(:));
  alpha   = repmat(alpha(:,i),1,ns);
  mu      = repmat(mu(:,i),1,ns);

  % Compute the unnormalized log-importance weights.
  fprintf('Computing importance weights for %d combinations ',ns);
  fprintf('of hyperparameters.\n');
  [logw alpha mu s] = outerloophyperbin2(X,y,A,B,alpha,mu,eta,h,log10odds);
