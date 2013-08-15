% SUMCOLS(X) returns a column vector in which each entry is the sum over
% each row of X. This is the same as SUM(X,2).
function y = sumcols (x)
  y = sum(x,2);