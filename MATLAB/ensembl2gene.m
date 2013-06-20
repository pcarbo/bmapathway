% ENSEMBL2GENE(FILE,IDS) creates a map from Ensembl access numbers to gene
% indices. FILE specifies the pathname to the 'gene2ensembl' file. IDS is a
% Map object (see CONTAINERS.MAP) specifying the mapping from Entrez gene ID
% to gene index. The return value is also a Map object.
function ensembl = ensembl2gene (file, ids)

  % Tab character.
  tab = char(9);  

  % Initialize a mapping from Ensembl accession numbers to Entrez gene
  % IDs. 
  ensembl = containers.Map('KeyType','double','ValueType','double');

  % Open the gene2ensembl file.
  fid = fopen(file,'r');
  
  % Skip the first line.
  str = fgetl(fid);

  % Repeat until we've reached the end of the file.
  while true
    
    % Read the next line from the file.
    str = fgetl(fid);
    if str == -1
      break
    end

    % Skip the first token (tax_id).
    [t str] = strtok(str,tab);

    % Get the Entrez gene ID.
    [t str] = strtok(str,tab);
    i = str2double(t);
  
    % Get the Ensembl accession number.
    t = strtok(str,tab);
    j = str2double(t(5:end));

    % If the gene ID is found, add a mapping from the Ensembl number to the
    % gene index.
    if isKey(ids,i)
      ensembl(j) = ids(i);
    end
  end

  % Close the file.
  fclose(fid);
