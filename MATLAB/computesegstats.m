function [P1, P2] = computesegstats (A, PIP, logw)

  % Get the number of SNPs (p) and the number of segments (n).
  [p n] = size(A);

  % Set up storage for the results.
  P1 = zeros(n,1);
  P2 = zeros(n,1);

  % Get the "weights" according to (composite) Simpson's rule.
  if isvector(logw)
    r = simpson(length(logw));
  else
    [n1 n2] = size(logw);
    r       = kron(simpson(n1),simpson(n2)');
  end

  % Repeat for each segment.
  for i = 1:n

    % Get the SNPs in the segment.
    snps = find(A(:,i));
      
    % Compute, for each setting of the hyperparameters, the posterior
    % probability that no SNPs in the segment are included in the model, and
    % the posterior probability that exactly one SNP in the segment is
    % included in the model.
    p0 = reshape(computeprob0(PIP(snps,:)),size(logw));
    p1 = reshape(computeprob1(PIP(snps,:)),size(logw));

    % Compute the posterior probability that at least 1 SNP in the segment is
    % included in the model (P1), and the posterior probability that at
    % least 2 SNPs in the segment are included in the model (P2), using
    % numerical integration to average over settings of the hyperparameters.
    P1(i) = quadpe(logw,1 - p0,r);
    P2(i) = quadpe(logw,1 - p1 - p0,r);
  end
