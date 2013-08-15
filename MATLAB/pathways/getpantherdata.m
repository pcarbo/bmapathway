% GETPANTHERDATA(FILE,ENSEMBL,N,M,NA) imports pathway data from a PANTHER
% pathway database file in "file in Gene Matrix Transposed format (.gmt), in
% which genes are specified by their Entrez gene ID. FILE is the pathname. M
% is the number of pathways to read (there is one pathway per line of the
% file). IDS is a Map object (see CONTAINERS.MAP) specifying the mapping
% from Entrez gene ID to gene index.
%
% The return value is a structure array with two fields, LABEL and GENES.
% LABEL is a cell array with M entries conntaining the pathway names. GENES
% is a sparse matrix, such that GENES(I,J) = 1 if and only if gene I is
% assigned to pathway J. Input N is the number of genes.
%
% NA specifies the maximum number of pathway annotations. If the number of
% pathway annotations (or equivalently the number of nonzero entries in
% GENES) is greater than NA, an error is reported.
function pathway = getpantherdata (file, ensembl, n, m, na)

  % Tab character.
  tab = char(9);  

  % Initialize the pathway information.
  pathway.label = cell(m,1);        % Pathway names.
  pathway.genes = spalloc(n,m,na);  % Annotation matrix.

  % Initialize the mapping from pathway accession to pathway index.
  pathids = containers.Map('KeyType','double','ValueType','double');  

  % This is the number of unique pathways encountered so far.
  t = 0;

  % This is the total number of gene-pathway annotations added so far to the
  % annotation matrix.
  a = 0;

  % Open the file.
  fid = fopen(file,'r');

  % Repeat until we've reached the end of the file.
  while true

    % Read the next line from the file.
    str = fgetl(fid);
    if str == -1
      break
    end

    % Get the pathway accession.
    [s str] = strtok(str,tab);
    pid = str2double(s(2:end));
    if isnan(pid)
      error('Not a pathway accession');
    end

    % Get the name of the pathway.
    [label str] = strtok(str,tab);

    % Check whether we have encountered this pathway.
    if ~isKey(pathids,pid)

      % This pathway has not been encountered in the database, so add it to
      % our collection. But first check whether we have reached the maximum
      % number of pathways.
      if t == m
	error('Reached the maximum number of pathways');
      end

      % Add the pathway to our pathway collection.
      t = t + 1;
      pathids(pid) = t;
      pathway.label{t} = label;
    end

    % Get the pathway index.
    j = pathids(pid);

    % Skip the next two tokens (pathway component accession, pathway
    % component name).
    [s str] = strtok(str,tab);
    [s str] = strtok(str,tab);

    % Get the Ensembl accession.
    s  = strtok(str,tab);
    k  = strfind(s,'ENSG');
    s  = s(k+4:end);
    s  = strtok(s,'|');
    id = str2double(s);
    if isnan(id)
      error('Not an Ensembl human gene accession');
    end

    % If this Ensembl accession number maps to a gene, add the annotation to
    % the annotation matrix.
    if isKey(ensembl,id)

      % Check whether we have reached the maximum number of annotations.
      if a == na
	error('Maximum number of annotations reached');
      end

      i = ensembl(id);
      pathway.genes(i,j) = 1;
      a = a + 1;
    end
  end

  % Close the file.
  fclose(fid);
