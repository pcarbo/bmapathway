% SUMROWS(X) returns a column vector in which each entry is the sum over
% each column of X. This is the same as TRANSPOSE(SUM(X,1)).
function y = sumrows (x)
  y = sum(x,1);