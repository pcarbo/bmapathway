% COMPUTEPROB0(PIP) returns the posterior probability that no variable is
% included in the model, assuming that the regression coefficients are
% independent of each other a posteriori. Each row of PIP corresponds to a
% variable, and each column corresponds to a different setting. The return
% value is a column vector equal to the number of columns of PIP.
function p0 = computeprob0 (PIP)
  p0 = prod(1 - PIP,1)';

  