% This is a script to load the pathway databases into MATLAB.
%
% NOTE: this script only works in MATLAB 7.10 (R2010a).
clear

% The number of pathways from each database.
numpc       = 1666;
numbiocarta = 301;
numpanther  = 140;
numbiosys   = 2235;

% The locations of the pathway database files.
pcfile         = 'homo-sapiens-9606-entrez-gene-id.gmt';
biocartafile   = 'biocarta.gmt';
pantherfile    = 'SequenceAssociationPathway3.01.txt';
biosyspathfile = 'biosystems.txt';
biosysgenefile = 'biosystems_gene';

% Load gene data.
fprintf('Loading gene ID, symbol and Ensembl accession maps.\n');
load('genemap.mat');

% Get the number of genes.
n = numvalues(ids);

% Load the Pathway Commons database.
fprintf('Importing Pathway Commons database.\n');
pc = getgmtdata(pcfile,ids,n,numpc,2e5);

% Load the BioCarta pathways.
fprintf('Importing BioCarta pathways.\n');
biocarta = getgmtdata2(biocartafile,symbols,n,numbiocarta,1e4);

% Load the PANTHER pathways.
fprintf('Importing PANTHER pathways.\n');
panther = getpantherdata(pantherfile,ensembl,n,numpanther,1e4);

% Load the BioSystems pathway database.
fprintf('Importing BioSystems pathway database.\n');
biosys = getbiosysdata(biosyspathfile,biosysgenefile,ids,n,numbiosys);

% Create two additional gene sets, one corresponding to the "classical"
% major histocompatibility complex (MHC), and another corresponding to the
% "extended" major histocompatibility complex (xMHC). Both are regions on
% the short arm of chromosome 6. The definitions of the MHC and xMHC are
% based on:
%
%   MHC Sequencing Consortium (1999) "Complete sequence and gene map of
%   a human major histocompatibility complex." Nature 401: 921–3.
%
%   Horton et al (2004) "Gene map of the extended human MHC." Nature 
%   Reviews Genetics 5: 889-899. 
%
% I define the gene sets for the MHC and xMHC to be the 124 and 251 gene
% loci, respectively, that are expressed within the MHC and xMHC regions.
%
mhcgenes  = cellfun(@(x) symbols(x),importdata('mhc_gene'));
xmhcgenes = cellfun(@(x) symbols(x),importdata('xmhc_gene'));
mhc.label = { 'Major histocompatibility complex'
                 'Extended major histocompatibility complex' };
mhc.genes = sparse(n,2);
mhc.genes(mhcgenes,1)  = 1;
mhc.genes(xmhcgenes,2) = 1;

% Merge pathways.
pc.database       = repmat({'PC'},numpc,1);
biocarta.database = repmat({'BioCarta'},numbiocarta,1);
biocarta.source   = repmat({'BioCarta'},numbiocarta,1);
panther.database  = repmat({'Panther'},numpanther,1);
panther.source    = repmat({'Panther'},numpanther,1);
biosys.database   = repmat({'BioSystems'},numbiosys,1);
mhc.database      = { 'NA'; 'NA' };
mhc.source        = { 'NA'; 'NA' };

pathway.label    = [ pc.label; biocarta.label; panther.label; 
                     biosys.label; mhc.label ];
pathway.database = [ pc.database; biocarta.database; panther.database;
		     biosys.database; mhc.database ];
pathway.source   = [ pc.source; biocarta.source; panther.source;
		     biosys.source; mhc.source ];
pathway.genes    = [ pc.genes biocarta.genes panther.genes ...
                     biosys.genes mhc.genes ];

% Remove pathways with no genes.
I = find(sum(pathway.genes) > 0);

pathway.label    = pathway.label(I);
pathway.database = pathway.database(I);
pathway.source   = pathway.source(I);
pathway.genes    = pathway.genes(:,I);

% Fix the source labels.
m = length(pathway.label);
for i = 1:m
  str = pathway.source{i};
  switch str
   case {'BIOCYC','HUMANCYC'}
    str = 'HumanCyc';
   case 'CELL_MAP'
    str = 'Cell Map';
   case {'NCI_NATURE','Pathway Interaction Database'}
    str = 'PID';
   case 'REACTOME'
    str = 'Reactome';
  end
  pathway.source{i} = str;
end

fprintf('Saving ungrouped pathway data.\n');
save('pathwayungrouped.mat','pathway','-v7.3');

% Group pathways into identical gene sets.
keep             = zeros(m,1);
pathway.synonyms = cell(m,1);
keep(1)          = 1;

% Repeat for each pathway except the first.
for i = 2:m
  
  % Check whether there is a pathway with the identical set of genes.
  I = 1:i-1;
  g = repmat(pathway.genes(:,i),1,length(I));
  j = find(keep(I) & sumrows(g ~= pathway.genes(:,I))' == 0);

  if length(j)

    % There is, so do not keep this pathway, and add it to the list of
    % "synonyms" for pathway j.
    pathway.synonyms{j} = [ pathway.synonyms{j}
		            struct('label',pathway.label{i},...
				   'database',pathway.database{i},...
				   'source',pathway.source{i}) ];
  else

    % There is not, so keep this pathway.
    keep(i) = 1;
  end
end

I = find(keep);

pathway.label    = pathway.label(I);
pathway.database = pathway.database(I);
pathway.source   = pathway.source(I);
pathway.synonyms = pathway.synonyms(I);
pathway.genes    = pathway.genes(:,I);

fprintf('Saving pathway data.\n');
save('pathway.mat','pathway','-v7.3');
