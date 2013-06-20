% PVE(SX,SA,LOG10ODDS) returns the prior estimate of the proportion of
% variance in a given outcome (e.g. a quantitative trait) that is explained
% by a collection of explanatory variables (e.g. genetic markers). SA is the
% prior variance of the additive effects on the outcome, LOG10ODDS is the
% (base 10) logarithm of the prior odds for inclusion. SX is the sum of the
% sample variances for all the explanatory variables. If inputs SA and
% LOG10ODDS are both not scalars, they must be numeric arrays of the same
% dimension.
function h = pve (sx, sa, log10odds)

  % This is the expected value of the sample genetic variance divided by
  % the variance of the residual.
  c = sa .* sigmoid10(log10odds) * sx;

  % This is the prior estimate of the proportion of variance in the outcome
  % that is explained by the variables. 
  h = c./(c + 1);
