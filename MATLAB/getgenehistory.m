% [IDS,SYMBOLS] = GETGENEHISTORY(FILE,IDS,SYMBOLS) creates a map from
% discontinued gene IDs to gene indices, and a map from symbols to gene
% indices. FILE specifies the pathname of the 'gene_history' file.
% 
% IDS and SYMBOLS are Map objects; see help for CONTAINERS.MAP. IDS is a
% mapping from gene identifiers to indices. SYMBOLS is a mapping from gene
% symbols to gene indices. On input, IDS must contain the mapping from
% current gene ID to gene index. For more information, see GETENTREZDATA.
function [ids, symbols] = getgenehistory (file, ids, symbols)

  % Tab character.
  tab = char(9);  
  
  % Open the gene history file for reading.
  fid = fopen(file,'r');
  
  % Skip the first line of the file.
  str = fgetl(fid);

  % Repeat for each discontinued gene ID (line in file).
  while true
    
    % Read the next line from the file.
    str = fgetl(fid);
    if str == -1
      break
    end

    % Skip the first token (tax_id).
    [t str] = strtok(str,tab);

    % Get the current gene identifier. If no identifier is provided, skip
    % this record.
    [t str] = strtok(str,tab);
    curid = str2double(t);
    if isfinite(curid)
  
      % Get the discontinued gene ID.
      [t str] = strtok(str,tab);
      oldid   = str2double(t);
    
      % Get the discontinued gene symbol.
      oldsymbol = strtok(str,tab);
    
      % Create a new entry in the gene map, provided that a match is found
      % in the existing gene map.
      if isKey(ids,curid)
	
	% Get the gene index.
	i = ids(curid);
	
	% Check whether the gene ID map is consistent. If not, report an
        % error. 
	if isKey(ids,oldid)
	  if ids(oldid) ~= i
	    error('Gene ID map is not consistent');
	  end
	end

	% Store the new mapping.
	ids(oldid) = i;

	% Check whether the gene symbol map is consistent. If so, there is
        % an ambiguous mapping from old symbol to new symbol. This is
        % unavoidable, so we retain the original mapping.
	if isKey(symbols,oldsymbol)
	  if symbols(oldsymbol) ~= i
	    i = symbols(oldsymbol);
	  end
	end
	
	% Store the new mapping.
	symbols(oldsymbol) = i;
      end
    end
  end
  
  % Close the file.
  fclose(fid);
