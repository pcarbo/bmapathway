% [GENE,IDS,SYMBOLS] = GETENTREZDATA(FILE,N) reads in gene information from
% an Entrez Gene 'gene_info' file, with pathname FILE. N is the number of
% genes (lines) to read in from the file.
%
% Output GENES is a structure array containing the following fields:
%
%    GENE.ID       gene identifiers
%    GENE.SYMBOL   gene symbols
%    GENE.CHR      chromosome numbers
%    GENE.DESC     description of gene
%
% Each field is an array with N entries. 
%
% Return values IDS and SYMBOLS are Map objects; for details, see help for
% CONTAINERS.MAP. IDS is a mapping from gene identifiers to indices, and
% SYMBOLS is a mapping from gene symbols to gene indices. For example, if
% GENE.ID(2105) = 2538 and GENE.SYMBOL(2105) = 'G6PC', then IDS(2538) and
% SYMBOLS('G6PC') both return 2105.
function [gene, ids, symbols] = getentrezdata (file, n)

  % If true, also load "unofficial" gene symbols.
  readunofficialsymbols = false;
  
  % Tab character.
  tab = char(9);  
  
  % Initialize the gene information.
  gene.id     = zeros(n,1);   % Gene identifier.
  gene.symbol = cell(n,1);    % Gene symbol.
  gene.chr    = zeros(n,1);   % Chromosome number.
  gene.desc   = cell(n,1);    % Gene description.

  % Initialize a map for the gene IDs and symbols.
  ids     = containers.Map('KeyType','double','ValueType','double');
  symbols = containers.Map('KeyType','char','ValueType','double');

  % Open the Entrez gene file for reading.
  fid = fopen(file,'r');
  
  % Repeat for each gene.
  for i = 1:n
  
    % Read the next line from the file.
    str = fgetl(fid);
    if str == -1
      error('Reached end of file');
    end

    % Skip the first token (tax_id).
    [t str] = strtok(str,tab);
  
    % Get the gene identifier.
    [t str] = strtok(str,tab);
    j = str2double(t);
    gene.id(i) = j;
  
    % Get the gene symbol. If the gene symbol has already been inserted
    % in the map, then there is an inconsistency, and report an error.
    [s str] = strtok(str,tab);
    if ~readunofficialsymbols & isKey(symbols,s)
      error('Gene symbols are inconsistent');
    end
    gene.symbol{i} = s;
    
    % Put the gene ID and symbol in the hash table.
    ids(j)     = i;
    symbols(s) = i;
    
    % Skip the next token (LocusTag).
    [t str] = strtok(str,tab);

    % Put the "unofficial" gene symbols into the hash table, but do not
    % override previously defined gene symbols.
    [t str] = strtok(str,tab);
    if readunofficialsymbols
      while length(t)
	[s t] = strtok(t,'|');
	if ~isKey(symbols,s)
	  symbols(s) = i;
	end
      end
    end
    
    % Skip the next token (dbXrefs).
    [t str] = strtok(str,tab);

    % Get the chromosome number.
    [t str] = strtok(str,tab);
    c = str2double(t);    
    if isfinite(c)
      gene.chr(i) = c;
    end
  
    % Skip the next token (map location).
    [t str] = strtok(str,tab);
  
    % Get the gene description.
    gene.desc{i} = strtok(str,tab);
  end

  % Close the file.
  fclose(fid);
  