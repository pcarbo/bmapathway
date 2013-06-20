% GETGMTDATA(FILE,IDS,N,M,NA) imports pathway data from a file in Gene
% Matrix Transposed format (.gmt), in which genes are specified by their
% Entrez gene ID. FILE is the pathname. M is the number of pathways to read
% (there is one pathway per line of the file). IDS is a Map object (see
% CONTAINERS.MAP) specifying the mapping from Entrez gene ID to gene index.
%
% The return value is a structure array with three fields, LABEL, SOURCE and
% GENES. LABEL and SOURCE are cell arrays with M entries. LABEL contains the
% pathway names, and SOURCE contains the pathway source (e.g. Reactome,
% Cancer Cell Map). GENES is a sparse matrix, such that GENES(I,J) = 1 if
% and only if gene I is assigned to pathway J. Input N is the number of
% genes. 
%
% NA specifies the maximum number of pathway annotations. If the number of
% pathway annotations (or equivalently the number of nonzero entries in
% GENES) is greater than NA, an error is reported.
function pathway = getgmtdata (file, ids, n, m, na)

  % Tab character.
  tab = char(9);  

  % Initialize the pathway information.
  pathway.label  = cell(m,1);        % Pathway names.
  pathway.source = cell(m,1);        % Pathway source.
  pathway.genes  = spalloc(n,m,na);  % Annotation matrix.

  % Open the file.
  fid = fopen(file,'r');

  % This is the total number of gene-pathway annotations added so far to the
  % annotation matrix.
  a = 0;

  % Repeat for each pathway.
  for j = 1:m

    % Read the next line from the pathway database file.
    str = fgetl(fid);
    if str == -1
      error('Reached end of file');
    end

    % Get the pathway label.
    [pathway.label{j} str] = strtok(str,tab);
    
    % Skip the pathway source.
    [pathway.source{j} str] = strtok(str,tab);  
  
    % The remaining tokens are the genes (identified by Entrez gene ID) that
    % participate in the pathway.
    while length(str)

      % Get the next token.
      [t str] = strtok(str,tab);  
    
      % Look up the gene by Entrez Gene ID.
      id = str2double(t);
      if isfinite(id)
	if isKey(ids,id)
	  i = ids(id);

	  % Check whether we have reached the maximum number of annotations.
	  if a == na
	    error('Maximum number of annotations reached');
	  end

	  % Add the annotation to the annotation matrix.
	  pathway.genes(i,j) = 1;
	  a = a + 1;
	end
      end
    end
  end

  % Close the file.
  fclose(fid);
