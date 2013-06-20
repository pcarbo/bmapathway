% Get the information about recombination rates and genetic distances
% across the genome from the HapMap data downloaded from the IMPUTE
% website. (The file is genetic_map_b35_combined.tgz.) 
clear

% Initialize the data.
chr  = [];  % Chromosome number.
pos  = [];  % Chromosomal position in base pairs.
freq = [];  % Recombination frequency (cM/Mb) across interval (i-1,i).
dist = [];  % Genetic distance (cM) from first position on chromsome.

% Repeat for each chromosome.
for c = 1:22
  fprintf('- Chromosome %d.\n',c);
  
  % Open the map file.
  data = importdata(sprintf('/tmp/pcarbo/genetic_map_chr%d.txt',c));
  data = data.data;
  
  % Get the number of genetic markers.
  n = size(data,1);
  
  % Get the chromosome numbers (chr), chromosomal positions (pos),
  % recombination frequencies (freq), and genetic distances (dist).
  chr  = [ chr;  repmat(c,n,1) ];
  pos  = [ pos;  data(:,1) ];
  freq = [ freq; data(:,2) ];
  dist = [ dist; data(:,3) ];
end

% Save the data.
hapmap = struct('chr',chr,'pos',pos,'freq',freq,'dist',dist);
save('/tmp/pcarbo/hapmap.mat','hapmap','-v7.3');
