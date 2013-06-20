clear

% Load the genotype and phenotype data.
load('/tmp/pcarbo/t1d.mat');

% Get the number of SNPs.
p = length(labels);

% OUTPUT CASES.
cases = find(y == 1);

% Open the file.
file = fopen('/tmp/pcarbo/t1d_cases.txt','w');

% Repeat for each SNP.
for j = 1:p
  fprintf('%6d',j);
  fprintf(repmat('\b',1,6));

  % Print the SNP label and the minor and major alleles.
  fprintf(file,'rs%d %c %c ',labels(j),minor(j),major(j));

  % Print the genotypes.
  fprintf(file,'%0.3f ',X(cases,j));
  fprintf(file,'\n');
end
fprintf('\n');

% Close the file.
fclose(file);

% OUTPUT CONTROLS.
controls = find(y == 0);

% Open the file.
file = fopen('/tmp/pcarbo/t1d_controls.txt','w');

% Repeat for each SNP.
for j = 1:p
  fprintf('%6d',j);
  fprintf(repmat('\b',1,6));

  % Print the SNP label and the minor and major alleles.
  fprintf(file,'rs%d %c %c ',labels(j),minor(j),major(j));

  % Print the genotypes.
  fprintf(file,'%0.3f ',X(controls,j));
  fprintf(file,'\n');
end
fprintf('\n');

% Close the file.
fclose(file);
