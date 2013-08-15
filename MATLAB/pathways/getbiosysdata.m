% GETBIOSYSDATA(PATHFILE,GENEFILE,IDS,N,M) imports pathway data from the
% NCBI BioSystems database. M is the number of pathways to read. IDS is a
% Map object (see CONTAINERS.MAP) specifying the mapping from Entrez gene ID
% to gene index.
%
% Two files are required to import the BioSystems pathway database. PATHFILE
% is the pathname of a file containing IDs and descriptors for the pathways
% in tab-delimited format. Each line of the file corresponds to a pathway.
% GENEFILE is the pathname of a file containing gene-pathway assignments.
% Each line of the file corresponds to an assignment of a gene to a pathway.
%
% The return value is a structure array with three fields, LABEL, SOURCE and
% GENES. LABEL and SOURCE are cell arrays with M entries. LABEL contains the
% pathway names, and SOURCE contains the pathway source (e.g. Reactome,
% KEGG). GENES is a sparse matrix, such that GENES(I,J) = 1 if and only if
% gene I is assigned to pathway J. Input N is the number of genes.
function pathway = getbiosysdata (pathfile, genefile, ids, n, m)

  % Tab character.
  tab = char(9);  

  % Initialize the pathway information.
  pathway.label  = cell(m,1);  % Pathway names.
  pathway.source = cell(m,1);  % Pathway source.

  % IMPORT PATHWAY NAMES AND IDS.
  % Initialize the mapping from pathway accession to pathway index.
  pathids = containers.Map('KeyType','double','ValueType','double');  

  % Open the pathway file.
  fid = fopen(pathfile,'r');
  
  % Skip the first line of the pathway file.
  str = fgetl(fid);
  
  % Repeat for each pathway.
  for j = 1:m
  
    % Read the next line from the pathway file.
    str = fgetl(fid);
    if str == -1
      error('Reached end of file');
    end
  
    % Get the pathway ID.
    [t str]      = strtok(str,tab);  
    pid          = str2double(t);
    pathids(pid) = j;
    
    % Skip the next token (accession).
    [t str] = strtok(str,tab);  
    
    % Get the pathway name.
    [pathway.label{j} str] = strtok(str,tab);  
    
    % Get the pathway source.
    [pathway.source{j} str] = strtok(str,tab);  
  end
  
  % Close the file.
  fclose(fid);
  
  % IMPORT PATHWAY-GENE ASSOCIATIONS.
  % Get the pathway and gene IDs.
  M   = importdata(genefile);
  pid = M(:,1);
  id  = M(:,2);
  pid = mat2cell(pid,ones(size(pid)),1);
  id  = mat2cell(id,ones(size(id)),1);
  clear M
  
  % Find the associations for which the pathway ID maps to a pathway, and
  % the Entrez gene ID maps to a gene. Note that some of the associations
  % may be redundant if two pathway IDs map to the same pathway, or if
  % two gene IDs map to the same gene.
  I   = find(isKey(ids,id) & isKey(pathids,pid));
  pid = pid(I);
  id  = id(I);

  % Create the annotation matrix.
  I = cell2mat(values(ids,id));
  J = cell2mat(values(pathids,pid));
  pathway.genes = spones(sparse(I,J,ones(size(I)),n,m));
