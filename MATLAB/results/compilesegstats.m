% For each segment, compute the posterior statistics P1 and P2 averaged over
% settings of the hyperparameters.
function [P1, P2] = compilesegstats (A, PIP, w)

  % Get the number of SNPs (p) and the number of segments (n).
  [p n] = size(A);

  % These vectors will contain the relevant posterior statistics for each
  % segment.
  P1 = zeros(n,1);
  P2 = zeros(n,1);

  % Repeat for each segment.
  for i = 1:n

    % Get the SNPs in the segment.
    snps = find(A(:,i));
      
    % Compute, for each setting of the hyperparameters, the posterior
    % probability that no SNPs in the segment are included in the model,
    % and the posterior probability that exactly one SNP in the segment is
    % included in the model.
    p0 = computeprob0(PIP(snps,:));
    p1 = computeprob1(PIP(snps,:));
      
    % Compute the posterior probability that at least 1 SNP in the segment is
    % included in the model (P1), and the posterior probability that at
    % least 2 SNPs in the segment are included in the model (P2), using a
    % simple numerical approximation to average over settings of the
    % hyperparameters.
    P1(i) = dot(w,1 - p0);
    P2(i) = dot(w,1 - p1 - p0);
  end
