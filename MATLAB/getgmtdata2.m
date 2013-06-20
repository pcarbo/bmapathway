% GETGMTDATA(FILE,SYMBOLS,N,M,NA) imports pathway data from a file in Gene
% Matrix Transposed format (.gmt), in which genes are specified by their
% gene symbol. FILE is the pathname. M is the number of pathways to read
% (there is one pathway per line of the file). SYMBOLS is a Map object (see
% CONTAINERS.MAP) specifying the mapping from gene symbol to gene index.
%
% The return value is a structure array with two fields, LABEL and GENES.
% LABEL is a cell array with M entries containing pathway names. GENES is a
% sparse matrix, such that GENES(I,J) = 1 if and only if gene I is assigned
% to pathway J. Input N is the number of genes.
%
% NA specifies the maximum number of pathway annotations. If the number of
% pathway annotations (or equivalently the number of nonzero entries in
% GENES) is greater than NA, an error is reported.
function pathway = getgmtdata2 (file, symbols, n, m, na)

  % Tab character.
  tab = char(9);  

  % Initialize the pathway information.
  pathway.label  = cell(m,1);        % Pathway names.
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

    % Skip the first token.
    [s str] = strtok(str,tab);  

    % Get the pathway label.
    [pathway.label{j} str] = strtok(str,tab);
  
    % The remaining tokens are the genes (identified by gene symbol) that
    % participate in the pathway.
    while length(str)

      % Get the next token.
      [s str] = strtok(str,tab);  

      % Look up the gene by gene symbol.
      if length(s)
	if isKey(symbols,s)
	  i = symbols(s);

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
