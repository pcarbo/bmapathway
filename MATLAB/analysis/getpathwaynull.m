function [logw, alpha, mu, s] = getpathwaynull (null, alt, theta0, theta, k)

  % Get the number of settings for the genome-wide log-odds (theta0) and
  % the enrichment parameter (n0).
  n0 = numel(alt.theta0);
  n1 = numel(alt.theta);

  % Round the settings of the genome-wide logodds to the nearest 0.01.
  theta0      = round(100*theta0)/100;
  null.theta0 = round(100*null.theta0)/100;
  alt.theta0  = round(100*alt.theta0)/100;

  % Get the statistics under the null hypothesis that match up with the
  % selected values of the genome-wide log-odds parameter (theta0).
  [ans I] = ismember(theta0,null.theta0);
  alpha   = null.alpha(:,I);
  mu      = null.mu(:,I);
  s       = null.s(:,I);

  % Get the indices corresponding to the selected settings of the
  % genome-wide log-odds (theta0) and the enrichment parameter (theta)
  % under the enrichment ("alternative") hypothesis.
  [ans I] = ismember(theta0,alt.theta0);
  j       = find(alt.theta == theta);

  % Get the log-importance weights corresponding to the selected settings of
  % the genome-wide log-odds, the chosen enrichment parameter setting, and
  % the chosen enrichment hypothesis.
  logw = alt.logw1(I,j,k);

  % Get the SNPs assigned to the enriched pathways.
  paths = find(alt.H(:,k));
  snps  = find(sumcols(alt.A(:,paths)));
  p     = length(snps);

  % Get the posterior statistics for these SNPs.
  alphapath = alt.alpha{k}(snps,:);
  mupath    = alt.mu{k}(snps,:);
  spath     = alt.s{k}(snps,:);
  alphapath = reshape(full(alphapath),[p n0 n1]);
  mupath    = reshape(full(mupath),[p n0 n1]);
  spath     = reshape(full(spath),[p n0 n1]);

  % Get the posterior statistics for all SNPs.
  alpha(snps,:) = alphapath(:,I,j);
  mu(snps,:)    = mupath(:,I,j);
  s(snps,:)     = spath(:,I,j);

