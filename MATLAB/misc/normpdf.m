% NORMPDF(X,MU,S) returns the probability density function of the univariate
% normal distribution at elements of X, with mean MU and variance S.
function p = normpdf (x, mu, s)
  p = exp(lognormal(x,mu,s));
