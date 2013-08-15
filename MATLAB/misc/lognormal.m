% LOGNORMAL(X,MU,S) returns the logarithm of the probability density
% function of the univariate normal distribution at elements of X, with mean
% MU and variance S.
function f = lognormal (x, mu, s)
  f = -log(2*pi*s)/2 - (x - mu).^2./(2*s);
