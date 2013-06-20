% [ANEW,CHRNEW,POSNEW] = COMBINESEGMENTS(A,CHR,POS) combines adjacent
% segments on the same chromosome, so that the newly combined segments
% overlap. For more information on the intput arguments, see the help for
% function MAKESEGMENTS.  The segments must be ordered by chromosome number,
% then by position along the chromosome.
function [Anew, chrnew, posnew] = combinesegments (A, chr, pos)

  % Get the number of segments.
  n = length(chr);

  % Get the number of chromosomes.
  c = max(chr);

  % Create the mapping from small segments to large segments.
  B      = spalloc(n,n-c,2*n);
  chrnew = zeros(n-c,1);
  posnew = zeros(n-c,1);

  % Repeat for each segment except the last. Here, i is the index of the
  % original segment, and j is the index of the combined segment.
  j = 0;
  for i = 1:n-1
    if chr(i) == chr(i+1)

      % Move to next the combined segment.
      j = j + 1;

      % Create the combined segment.
      B(i,j)    = 1;
      B(i+1,j)  = 1;
      chrnew(j) = chr(i);
      posnew(j) = pos(i);
    end
  end

  % Create the new adjacency matrix.
  Anew = A * B;
