% VAROFMIXTURE(MU,S,W) computes the variance of random vector X, in which X
% is a mixture of multivariate normals with means MU and variances S. Inputs
% MU and S are matrices of dimension N x K, where N is the dimension of the
% random variable X, and K is the number of components in the mixture;
% MU(:,i) and S(:,i) are the mean and variance of mixture component i.
% Input W is an array containing the mixture weights W(i) for each mixture
% component i. These mixture weights need not be normalized.
function y = varofmixture (mu, s, w)

  % Normalize the mixture weights.
  w = w(:);
  w = w / sum(w);

  % Compute the second moment E(X^2) for each mixture component.
  E2 = s + mu.^2;

  % Compute the variance of X.
  y = E2*w - (mu*w).^2;
