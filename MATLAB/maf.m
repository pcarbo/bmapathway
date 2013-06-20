% MAF(X) returns maximum likelihood estimates of minor allele frequencies
% for genotype data X. Input X is an N x P matrix of minor allele counts,
% where N is the number of individuals (or samples), and P is the number of
% genetic loci (or SNPs).
function f = maf (X)
  f = mean(X)'/2;

  