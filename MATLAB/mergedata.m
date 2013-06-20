% [X,Y,LABELS,CHR,POS,MINOR,MAJOR] = MERGEDATA(X1,Y1,LABELS1,CHR1,...
% POS1,MINOR1,MAJOR1,X2,Y2,LABELS2,CHR2,POS2,MINOR2,MAJOR2) merges genotype
% and phenotype data from two data sets. It only merges data for SNPs that
% are the same in both data sets, according to refSNP identifiers LABELS1
% and LABELS2. After merging, it adjusts the alleles so that all minor
% allele frequencies are less than 0.5. When merging the samples, it also
% randomizes the order of the samples.
function [X, y, labels, chr, pos, minor, major] = ...
      mergedata (X1, y1, labels1, chr1, pos1, minor1, major1, ...
		 X2, y2, labels2, chr2, pos2, minor2, major2)

  % Get the number of individuals (samples) in each data set.
  n1 = length(y1);
  n2 = length(y2);
  n  = n1 + n2;

  % Find the SNPs shared by both data sets.
  [ans I J] = intersect(labels1,labels2);

  % Sort the SNPs from the first data set.
  X1      = X1(:,I);
  labels1 = labels1(I);
  chr1    = chr1(I);
  pos1    = pos1(I);
  minor1  = minor1(I);
  major1  = major1(I);

  % Sort the SNPs from the first data set.
  X2      = X2(:,J);
  labels2 = labels2(J);
  chr2    = chr2(J);
  pos2    = pos2(J);
  minor2  = minor2(J);
  major2  = major2(J);

  % Get the refSNP identifiers, the chromsome numbers, the
  % chromosomal positions, and the major and minor alleles.
  labels  = labels1;
  chr     = chr1;
  pos     = pos1;
  minor   = minor1;
  major   = major1;

  % Initialize storage for the merged genotypes. I store the genotypes in
  % single precision to use less memory.
  p = length(labels);
  X = zeros(n,p,'single');

  % Find SNPs for which minor alleles match.
  A = (minor1 == minor2 & minor1 ~= 'N' & minor2 ~= 'N') | ...
      (major1 == major2 & major1 ~= 'N' & major2 ~= 'N');
  I = find(A);

  % Merge SNPs for which minor alleles match.
  X(1:n1,I)   = X1(:,I);
  X(n1+1:n,I) = X2(:,I);

  % Find SNPs for which major alleles match minor alleles.
  B = (minor1 == major2 & minor1 ~= 'N' & major2 ~= 'N') | ...
      (major1 == minor2 & major1 ~= 'N' & minor2 ~= 'N');
  I = find(B);

  % Merge SNPs for which major alleles match minor alleles. Here we need to
  % the set the genotypes X for the second set of subjects to 2 - X.
  X(1:n1,I)   = X1(:,I);
  X(n1+1:n,I) = 2 - X2(:,I);
  clear X1 X2
  
  % If there are any SNPs for which neither minor or major alleles match,
  % report an error.
  if sum(~(A | B))
    error('Some SNPs do not have matching alleles');
  end

  % Make sure the minor alleles all have frequency less than 0.5.
  I        = find(maf(X) > 0.5);
  X(:,I)   = 2 - X(:,I);
  minor(I) = major1(I);
  major(I) = minor1(I);

  % Sort the SNPs by position.
  [ans I] = sort(chr * 1e9 + pos);
  X       = X(:,I);
  labels  = labels(I);
  pos     = pos(I);
  chr     = chr(I);
  minor   = minor(I);
  major   = major(I);

  % Get the phenotypes.
  y = [ y1; y2 ];
  
  % Randomize the order of the samples.
  is = randperm(n);
  X  = X(is,:);
  y  = y(is);
