% [A,SEGCHR,SEGPOS] = MAKESEGMENTS(CHR,POS,D,N) creates contiguous genomic
% segments for genetic loci. Inputs CHR and POS are vectors of length P;
% these vectors specify the chromosome number and position along the
% chromosome (in bases) for each SNP. It is assumed that the SNPs are
% ordered by chromosome number, then by position along the chromosome. D is
% the length of each segment (in bases), and N is the maximum number of
% segments to create.
%
% Output A is a sparse matrix of dimension P x M, where P is the number of
% SNPs and M is the number of segments created. A(I,J) = 1 if and only if
% genetic locus I belongs to the segment J. SEGCHR and SEGPOS are the
% chromosome number and start position of each segment.
function [A, segchr, segpos] = makesegments (chr, pos, d, n)

  % Get the number of SNPs.
  p = length(chr);

  % Initialize the outputs.
  A      = spalloc(p,n,p);  % Adjacency matrix.
  segchr = zeros(n,1);      % Chromosome number. 
  segpos = zeros(n,1);      % Start position of segment.

  % Repeat for each SNP. Here, i indexes segments and j indexes SNPs.
  i = 1;
  segpos(1) = pos(1);
  for j = 1:p

    % If the position of the current SNP is more than distance D away from
    % the start of the current segment, or if the current SNP is in a new
    % chromosome, create a new segment.
    if pos(j) > segpos(i) + d | chr(j) > segchr(i)

      % If we have already reached the maximum number of segments, stop.
      if i == n
	break
      end

      % Create a new segment.
      i = i + 1;
      segchr(i) = chr(j);
      segpos(i) = pos(j);
    end
    
    % Include this SNP in the current segment.
    A(j,i) = 1;
  end

  % Remove segments with no SNPs.
  I      = find(full(sumrows(A)));
  A      = A(:,I);
  segchr = segchr(I);
  segpos = segpos(I);
