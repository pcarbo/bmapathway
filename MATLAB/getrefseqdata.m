% [START,STOP] = GETREFSEQDATA(FILE,IDS,N) reads in gene sequence
% information from a RefSeq assembly ('seq_gene') file. FILE specifies the
% pathname of the RefSeq assembly file. N is the number of genes. Outputs
% START and STOP are the location of the 5' end and 3' end, respectively,
% for each gene. If no position information is available for a gene
% (i.e. if a gene does not map to RefSeq, the reference Sequence), then
% the respective entries of START and STOP are set to -1.
% 
% Input IDS must be be a mapping from current and discontinued gene
% identifiers to gene indices. For more information on this, see
% GETENTREZDATA and GETGENEHISTORY. This is needed in case the RefSeq
% assembly file contains outdated gene IDs.
function [start, stop] = getrefseqdata (file, ids, n)
  
  % Tab character.
  tab = char(9);   

  % Initialize the gene information. START is the location of 5' end, and
  % STOP is the location o the 3' end.
  start = repmat(-1,n,1);
  stop  = repmat(-1,n,1);
  
  % Open the RefSeq gene file for reading.
  fid = fopen(file,'r');
  
  % Skip the first line of the file.
  str = fgetl(fid);

  % Repeat for each line in the RefSeq gene file.
  while true

    % Read the next line from the file.
    str = fgetl(fid);
    if str == -1
      break
    end

    % Skip the first two tokens (taxid and chromosome).
    [t str] = strtok(str,tab);
    [t str] = strtok(str,tab);

    % Get the gene's start and stop positions on the chromosome.
    [t str] = strtok(str,tab);
    p0      = str2double(t);
    [t str] = strtok(str,tab);
    p1      = str2double(t);
  
    % Skip the next 6 columns (chr_orient, contig, ctg_start, ctg_stop,
    % ctg_orient and feature_name).
    [t str] = strtok(str,tab);
    [t str] = strtok(str,tab);
    [t str] = strtok(str,tab);
    [t str] = strtok(str,tab);
    [t str] = strtok(str,tab);
    [t str] = strtok(str,tab);

    % Get the gene ID.
    t  = strtok(str,tab);
    id = str2double(t(8:end));
    
    % Find the gene index.
    if isKey(ids,id)

      % An entry was found, so store the sequence information.
      i        = ids(id);
      start(i) = p0;
      stop(i)  = p1;
    end
  end
  
  % Close the file.
  fclose(fid);
  