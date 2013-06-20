% SPDIAG(X) returns a sparse diagonal matrix, with diagonal entries formed
% by vector X.
function X = spdiag (x)
  X = diag(sparse(x));
  