% [X,LABELS,CHR,POS,MINOR,MAJOR] = READBIMBAM2(N,P,GENFILES,POSFILES)
% imports the genotype data in BIMBAM format, in which the genotype data for
% each SNP is given by two numbers, the probability of minor allele count
% zero, and the probability of minor allele count one.
%
% Input N specifies the number of samples, and P specifies the number of
% SNPs. Input argument GENFILES is a string specifying the pathnames of the
% BIMBAM files containing the genotype data. Input POSFILES specifies the
% location of the files containing the chromosomal position data for the
% SNPs. The strings GENFILES and POSFILES must contain a placeholder
% (e.g. '%d') for the chromosome number.
%
% The six outputs are:
%
%     X       N x P matrix of mean genotypes
%     LABELS  vector of SNP identifiers
%     CHR     chromosome numbers
%     POS     chromosomal positions
%     MINOR   minor alleles
%     MAJOR   major alleles
%
% The position and genotype files must contain information for a total of P
% SNPs. The SNPs must be listed in the same order in the genotype and
% position files.
function [X, labels, chr, pos, minor, major] = ...
      readbimbam2 (n, p, genfiles, posfiles)

  % The number of chromosomes.
  N = 22;

  % Allocate storage for the data. I store the genotypes in single
  % precision to use less memory.
  X      = zeros(n,p,'single');
  labels = zeros(p,1);       % rs number from position files.
  chr    = zeros(p,1);       % Chromosome number.
  pos    = zeros(p,1);       % Chromosomal position.
  minor  = repmat(' ',p,1);  % Minor allele.
  major  = repmat(' ',p,1);  % Major allele.

  % This is the SNP index.
  ks = 0;

  % Repeat for each chromosome.
  fprintf('Reading SNP data.\n');
  for c = 1:N
    fprintf('Chromosome %02d',c);
    fprintf(repmat('\b',1,13));

    % Import the position and genotype data.
    clear a d
    a = importdata(sprintf(posfiles,c));
    d = importdata(sprintf(genfiles,c));
    
    % Get the number of SNPs on the chromosome.
    m  = length(a.data);
    ks = ks(end) + (1:m);

    % Get the chromosomal position and chromosome number.
    chr(ks) = c;
    pos(ks) = str2double(a.textdata(:,2));
    
    % Get the major and minor alleles.
    minor(ks) = cell2mat(d.textdata(:,2));
    major(ks) = cell2mat(d.textdata(:,3));        

    % Get genotype probabilities: P0 is the probability that the minor
    % allele count is 0; P1 is the probability that the minor allele is 1;
    % and P2 is the probability that the minor allele count is 2.
    p0 = d.data(:,1:2:end)';
    p1 = d.data(:,2:2:end)';
    p2 = 1 - p1 - p0;

    % The mean genotype (minor allele count) is simply P1 + 2*P2.
    X(:,ks) = p1 + 2*p2;

    % Get the SNP identifiers.
    labels(ks) = cellfun(@(x) getsnplabel(x),a.textdata(:,1));
  end
  fprintf('\n');
  
% ------------------------------------------------------------------
% Convert the label from a string to an integer. If the label is positive,
% then this is the refSNP identifier.
function x = getsnplabel (s)
  if strcmp('rs',s(1:2))
    x = str2double(s(3:end));
  else
    x = 0;
  end
