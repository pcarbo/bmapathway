% This is a script to load the relevant gene information into MATLAB for
% human genome assembly reference versions 17 (NCBI build 35.1) and 18 (NCBI
% build 36.3).
clear

% Number of genes.
n = 42097;   

% Load the information from the Entrez gene file.
fprintf('Reading Entrez gene data.\n');
[gene ids symbols] = getentrezdata('Homo_sapiens.gene_info',n);

% There should be at least one symbol and one ID for each gene. If not,
% report an error.
if numvalues(ids) < n
  error('Some genes do not map from Entrez gene IDs');
end
if numvalues(symbols) < n
  error('Some genes do not map from gene symbols');
end

% Load discontinued gene IDs and symbols.
fprintf('Reading gene history data.\n');
[ids symbols] = getgenehistory('gene_history',ids,symbols);

% Load mapping for Ensembl accession numbers.
fprintf('Reading gene2ensembl data.\n');
ensembl = ensembl2gene('gene2ensembl',ids);

% Save the gene identifier, symbol and Ensembl ID maps.
fprintf('Saving gene ID, symbol and Ensembl ID maps.\n');
save('genemap.mat','ids','symbols','ensembl','-v7.3');

% Load the gene positions for Build 35.1.
fprintf('Reading RefSeq data for build 35.1.\n');
[gene.start gene.stop] = getrefseqdata('seq_gene_35.1.md',ids,n);

% Save the gene data for Build 35.1.
fprintf('Saving gene data for build 35.1.\n');
save('gene_35.1.mat','gene','-v7.3');

% Load the gene positions for Build 36.3.
fprintf('Reading RefSeq data for build 36.3.\n');
[gene.start gene.stop] = getrefseqdata('seq_gene_36.3.md',ids,n);

% Save the gene data for Build 36.3.
fprintf('Saving gene data for build 36.3.\n');
save('gene_36.3.mat','gene','-v7.3');

