% GETNULLSTATS(FILE,THETA0) loads posterior statistics from FILE, and
% returns a STRUCT containing the posterior statistics under the null
% hypothesis that match up with the selected values or the genome-wide
% log-odds parameter (THETA0).
function stats = getnullstats (file, theta0)

  % Load the posterior statistics from file.
  stats = load(file);

  % Round the settings of the genome-wide logodds to the nearest 0.01.
  theta0       = round(100*theta0)/100;
  stats.theta0 = round(100*stats.theta0)/100;

  % Select the posterior statistics that match up with the selected
  % values for the genome-wide log-odds.
  [ans I]      = ismember(theta0,stats.theta0);
  stats.theta0 = stats.theta0(I);
  stats.logw   = stats.logw(I);
  stats.alpha  = stats.alpha(:,I);
  stats.mu     = stats.mu(:,I);
  stats.s      = stats.s(:,I);
  stats.eta    = stats.eta(:,I);
 

