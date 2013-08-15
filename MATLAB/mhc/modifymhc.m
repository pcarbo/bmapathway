% MODIFYMHC(PATHWAY,A) returns modified pathway annotations such that the
% annotations for the major histocompatibiltiy complex (MHC) and extended
% MHC include all SNPs in the region, not just the ones near the (known)
% expressed genes.
function A = modifymhc (pathway, A)

  % Modify the SNP annotations corresponding to these gene sets.
  labels = { 'Major histocompatibility complex'
             'Extended major histocompatibility complex' };

  n = length(labels);
  for i = 1:n

    % Get the SNPs assigned to the gene set.
    i    = find(strcmpi(pathway.label,labels{i}));
    snps = find(A(:,i));

    % Modify the annotation so that it includes all SNPs between the
    % first and last SNPs.
    a        = min(snps);
    b        = max(snps);
    A(a:b,i) = 1;
  end
