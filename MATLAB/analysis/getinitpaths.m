% GETINITPATHS(GENE,PATHWAY) returns a sparse matrix that specifies the
% initial set of enrichment hypotheses. I ignore pathways with less than one
% gene assigned to the reference genome. I also ignore NCI/Nature PID
% pathways from Pathway Commons with more than 500 genes.
function H = getinitpaths (gene, pathway)

  % Get the total number of candidate pathways or gene sets.
  m = length(pathway.label);

  % Get the number of genes assigned to each gene set.
  numgenes = pathway.genes' * inref(gene);

  % Remove pathways with less than one gene assigned to the reference
  % genome, and remove NCI/Nature PID pathways from Pathway Commons with
  % more than 500 genes.
  bad   = strcmp(pathway.database,'PC') & ...
          strcmp(pathway.source,'PID') & ...
          numgenes > 500;
  paths = find(numgenes > 1 & ~bad);

  % Sort the pathways by the number of genes in their gene set so that the
  % Bayes factors are computed for the smallest pathways first.
  [ans I] = sort(numgenes(paths));
  paths   = paths(I);
  H       = speye(m);
  H       = H(:,paths);
