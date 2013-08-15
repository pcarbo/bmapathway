% SNPS2GENES(CHR,POS,GENE,D,M) assigns SNPs to genes. Inputs CHR and POS
% given the chromosome number and position along the chromosome (in bases)
% for each SNP. Input GENE is a structure array with at least these three
% fields: GENE.CHR, the chromosome number for each gene; GENE.START, the
% location of the 5' end, in bases; and GENE.STOP, the location of the 3'
% end. Input D is the length of the "window" for assigning SNPs to genes: a
% SNP is assigned to a gene if it is within D bases of the gene's
% transcribed region. 
%
% Output A is a sparse matrix of assignments, such that A(I,J) = 1 if and
% only if SNP I is assigned to gene J. Input M specifies the maximum number
% of assignments, which is the same as the number of nonzeros in the
% assignment matrix. If SNPS2GENES attempts to make more than M assignments,
% an error is reported.
function A = snps2genes (chr, pos, gene, d, m)

  % Get the number of SNPs (p) and the number of genes (n).
  p = length(chr);
  n = length(gene.chr);

  % Initialize the sparse assignment matrix.
  A = spalloc(p,n,m);

  % This is the total number of assignments so far.
  a = 0;

  % Repeat for each gene.
  for i = 1:n

    % Check whether a chromosomal position has been assigned to the gene.
    if gene.start(i) >= 0 & gene.stop(i) >= 0

      % Find all SNPs that are located within D bases of the transcribed
      % region of the gene.
      snps = find(chr == gene.chr(i) & ...
		  pos > gene.start(i) - d & ...
		  pos < gene.stop(i) + d);

      % Check whether we have reached the maximum number of assignments.
      if a + length(snps) > m
	error('Maximum number of assignments reached');
      end

      % Add the SNP-gene assignments to the matrix.
      A(snps,i) = 1;
      a = a + length(snps);
    end
  end

