% DIAGPROD(A,B) efficiently computes DIAG(A*B).
function y = diagprod (A, B)
  y = sumcols(A.*B');
