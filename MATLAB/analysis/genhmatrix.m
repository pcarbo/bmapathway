% GENHMATRIX(PATHS,M) returns an M x N sparse matrix, where M is the number
% of candidate pathways (or gene sets), and N = NUMEL(PATHS) is the number
% of enrichment hypotheses. Input PATHS is a cell array in which each entry
% specifies the indices of the enriched pathways in the corresponding
% enrichment hypothesis. The return value H has nonzero entries H(j,i) = 1
% if and only if pathway j is enriched for associations with the binary
% trait in enrichment hypothesis i; otherwise, H(j,i) = 0.
function H = genhmatrix (paths, m)

  % Get the number of enrichment hypotheses.
  n = numel(paths);

  % Allocate space for the sparse matrix.
  H = spalloc(m,n,sum(cellfun(@(x)numel(x),paths)));

  % Add the enrichment hypotheses to the matrix.
  for i = 1:n
    H(paths{i},i) = 1;
  end

