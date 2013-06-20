% COMPUTEPROB1(PIP) returns the posterior probability that exactly one SNP
% is included in the model, assuming that the regression coefficients are
% independent of each other a posteriori. Each row of PIP corresponds to a
% variable, and each column corresponds to a different setting. The return
% value is a column vector equal to the number of columns of PIP.
function p1 = computeprob1 (PIP)
  p1 = computeprob0(PIP) .* sumrows(PIP./(1-PIP+eps))';
