% This is a script to import the WTCCC controls, bipolar disorder (BD)
% cases, coronary artery disease (CAD) cases, Crohn's disease (CD) cases,
% hypertension (HT) cases, rheumatoid arthritis (RA) cases, type 1 diabetes
% (T1D) cases, and type 2 diabetes (T2D) cases from files in BIMBAM format,
% in which the genotype data for each chromosome is stored in a separate
% file.

% 1958 BIRTH CONTROLS (58BC)
% --------------------------
clear

n = 1480;    % Number of subjects.
p = 459446;  % Number of SNPs.

posfiles = 'birthctrl.wtccc.chr%d.pos';
genfiles = 'birthctrl.chr%d.mean.genotype.txt';
matfile  = 'wtccc_58bc.mat';

fprintf('1958 birth controls\n');
[X labels chr pos minor major] = readbimbam(n,p,genfiles,posfiles);
save(matfile,'X','labels','chr','pos','minor','major','-v7.3');

% UK BLOOD SERVICES CONTROLS (UKBS)
% ---------------------------------
clear

n = 1458;    % Number of subjects.
p = 459446;  % Number of SNPs.
  
posfiles = 'bloodctrl.wtccc.chr%d.pos';
genfiles = 'bloodctrl.chr%d.mean.genotype.txt';
matfile  = 'wtccc_ukbs.mat';

fprintf('UK blood services controls\n');
[X labels chr pos minor major] = readbimbam(n,p,genfiles,posfiles);
save(matfile,'X','labels','chr','pos','minor','major','-v7.3');

% BIPOLAR DISORDER CASES (CD)
% ---------------------------
clear

n = 1868;    % Number of subjects.
p = 458868;  % Number of SNPs.

posfiles = 'bd.wtccc.chr%d.pos';
genfiles = 'bd.chr%d.mean.genotype.txt';
matfile  = 'wtccc_bd.mat';

fprintf('Bipolar disorder cases\n');
[X labels chr pos minor major] = readbimbam(n,p,genfiles,posfiles);
save(matfile,'X','labels','chr','pos','minor','major','-v7.3');

% CORONARY ARTERY DISEASE CASES (CAD)
% -----------------------------------
clear

n = 1926;    % Number of subjects.
p = 458868;  % Number of SNPs.

posfiles = 'cad.wtccc.chr%d.pos';
genfiles = 'cad.chr%d.mean.genotype.txt';
matfile  = 'wtccc_cad.mat';

fprintf('Coronary artery disease cases\n');
[X labels chr pos minor major] = readbimbam(n,p,genfiles,posfiles);
save(matfile,'X','labels','chr','pos','minor','major','-v7.3');

% CROHN'S DISEASE CASES (CD)
% --------------------------
clear

n = 1748;    % Number of subjects.
p = 458868;  % Number of SNPs.

posfiles = 'cd.wtccc.chr%d.pos';
genfiles = 'cd.chr%d.mean.genotype.txt';
matfile  = 'wtccc_cd.mat';

fprintf('Crohn''s disease cases\n');
[X labels chr pos minor major] = readbimbam(n,p,genfiles,posfiles);
save(matfile,'X','labels','chr','pos','minor','major','-v7.3');

% HYPERTENSION CASES (HT)
% -----------------------
clear

n = 1952;    % Number of subjects.
p = 458868;  % Number of SNPs.

posfiles = 'ht.wtccc.chr%d.pos';
genfiles = 'ht.chr%d.mean.genotype.txt';
matfile  = 'wtccc_ht.mat';

fprintf('Hypertension cases\n');
[X labels chr pos minor major] = readbimbam(n,p,genfiles,posfiles);
save(matfile,'X','labels','chr','pos','minor','major','-v7.3');

% RHEUMATOID ARTHRITIS (RA)
% -------------------------
clear

n = 1860;    % Number of subjects.
p = 458868;  % Number of SNPs.

posfiles = 'ra.wtccc.chr%d.pos';
genfiles = 'ra.chr%d.mean.genotype.txt';
matfile  = 'wtccc_ra.mat';

fprintf('Rheumatoid arthritis cases\n');
[X labels chr pos minor major] = readbimbam(n,p,genfiles,posfiles);
save(matfile,'X','labels','chr','pos','minor','major','-v7.3');

% TYPE 1 DIABETES CASES (T1D)
% ---------------------------
clear

n = 1963;    % Number of subjects.
p = 458868;  % Number of SNPs.
  
posfiles = 't1d.wtccc.chr%d.pos';
genfiles = 't1d.chr%d.mean.genotype.txt'; 
matfile  = 'wtccc_t1d.mat';

fprintf('Type 1 diabetes cases\n');
[X labels chr pos minor major] = readbimbam(n,p,genfiles,posfiles);
save(matfile,'X','labels','chr','pos','minor','major','-v7.3');

% TYPE 2 DIABETES (T2D)
% ---------------------
clear

n = 1924;    % Number of subjects.
p = 458868;  % Number of SNPs.

posfiles = 't2d.wtccc.chr%d.pos';
genfiles = 't2d.chr%d.mean.genotype.txt';
matfile  = 'wtccc_t2d.mat';

fprintf('Type 2 diabetes cases\n');
[X labels chr pos minor major] = readbimbam(n,p,genfiles,posfiles);
save(matfile,'X','labels','chr','pos','minor','major','-v7.3');
