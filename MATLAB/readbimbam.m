% [X,LABELS,CHR,POS,MINOR,MAJOR] = READBIMBAM(N,P,GENFILES,POSFILES) imports
% the genotype data in BIMBAM format. Input N specifies the number of
% samples, and P specifies the number of SNPs. Input argument GENFILES is a
% string specifying the pathnames of the BIMBAM files containing the
% genotype data. Input POSFILES specifies the location of the files
% containing the chromosomal position data for the SNPs. The strings
% GENFILES and POSFILES must contain a placeholder (e.g. '%d') for the
% chromosome number.
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
% These outputs are sorted by chromosome number, then by position along the
% chromosome. 
%
% The position and genotype files must contain information for a total of P
% SNPs. The SNPs may be listed in different orders in the genotype and
% position files.
function [X, labels, chr, pos, minor, major] = ...
      readbimbam (n, p, genfiles, posfiles)

  % The number of chromosomes.
  N = 22;

  % Allocate storage for the data. I store the genotypes in single
  % precision to use less memory.
  X       = zeros(n,p,'single');
  labels  = zeros(p,1);       % rs number from position files.
  labels2 = zeros(p,1);       % rs number from the genotype files.
  chr     = zeros(p,1);       % Chromosome number.
  pos     = zeros(p,1);       % Chromosomal position.
  minor   = repmat(' ',p,1);  % Minor allele.
  major   = repmat(' ',p,1);  % Major allele.

  % This is the SNP index.
  snps = 0;
  
  % Repeat for each chromosome.
  fprintf('Reading SNP data.\n');
  for c = 1:N
    fprintf('Chromosome %02d',c);
    fprintf(repmat('\b',1,13));

    % IMPORT POSITION DATA.
    filename = sprintf(posfiles,c);
    fid      = fopen(filename,'r');
    a        = importdata(sprintf(posfiles,c),'\t');

    % Get the number of SNPs on the chromosome.
    m    = length(a.textdata);
    snps = snps(end) + (1:m);

    % Get the refSNP identifier, chromosome number, and position on the
    % chromosome. 
    labels(snps) = cellfun(@(x) getsnplabel(x),a.textdata);
    pos(snps)    = a.data(:,1);
    chr(snps)    = a.data(:,2);

    % IMPORT GENOTYPE DATA.
    fid = fopen(sprintf(genfiles,c),'r');
    a   = textscan(fid,['%s %c %c ',repmat('%f ',1,n)]);
    fclose(fid);

    % Get the refSNP identifier, the major and minor alleles, and the mean
    % genotypes. The second and third columns of the genotype file are the
    % minor and major alleles, respectively.
    labels2(snps) = cellfun(@(x) getsnplabel(x),a{1});
    minor(snps)   = a{2};
    major(snps)   = a{3};
    X(:,snps)     = cell2mat(a(4:end))';
  end
  fprintf('\n');

  % Sort the chromosome position data by rs number.
  fprintf('Sorting SNPs.\n');
  [ans I] = sort(labels);
  labels  = labels(I);
  chr     = chr(I);
  pos     = pos(I);

  % Sort the genotype data by rs number.
  [ans I] = sort(labels2);
  labels2 = labels2(I);
  minor   = minor(I);
  major   = major(I);
  X       = X(:,I);

  % Sort everything by chromsome number, then by chromsomal position.
  [ans I] = sort(chr*1e9 + pos);
  labels  = labels(I);
  pos     = pos(I);
  chr     = chr(I);
  minor   = minor(I);
  major   = major(I);
  X       = X(:,I);
  
% ------------------------------------------------------------------
% Convert the label from a string to an integer. If the label is positive,
% then this is the refSNP identifier.
function x = getsnplabel (s)
  if strcmp('rs',s(1:2))
    x = str2double(s(3:end));
  elseif strcmp('SNP_A-',s(1:6))
    x = -str2double(s(7:end));
  else
    x = 0;
  end
  