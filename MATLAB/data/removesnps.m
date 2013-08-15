% [X,LABELS,CHR,POS,MINOR,MAJOR] = REMOVESNPS(X,LABELS,CHR,POS,...
% MINOR,MAJOR,SNPS) removes the indices specified by SNPS from the data set,
% and retains the others.
function [X, labels, chr, pos, minor, major] = ...
      removesnps (X, labels, chr, pos, minor, major, snps)
  
  % Get the set of SNPs to keep.
  p    = length(labels);
  snps = setdiff(1:p,snps);
  
  % Keep the specified SNPs.
  X      = X(:,snps);
  labels = labels(snps);
  chr    = chr(snps);
  pos    = pos(snps);
  minor  = minor(snps);
  major  = major(snps);
  