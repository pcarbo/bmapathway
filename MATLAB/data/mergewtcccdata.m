% Creates the final case-control data sets for bipolar disorder (BD),
% coronary artery disease (CAD), Crohn's disease (CD), hypertension (HT),
% rheumatoid arthritis (RA), type 1 diabetes (T1D) and type 2 diabetes
% (T2D). Run this script after running the GETWTCCCDATA script on the
% Wellcome Trust case-control data.
%
% For all seven diseases, I have an additional quality control step to
% filter out potentially problematic SNPs. Some common SNPs (i.e. that have
% high minor allele frequencies) show moderate evidence for association
% based on my calculations (single-SNP Bayes factors with the prior standard
% deviation of the log-odds ratios set to 0.1), but they do not match up
% with any of the SNPs reported in Supplementary Table 7 of the original
% publication (this table lists SNPs with trend p-values less than 1e-4 that
% have at least one other nearby SNP with a p-value less than 1e-3, and
% because they do not appear to be supported by nearby SNPs based on
% inspecting the single-SNP Bayes factors, we cannot rule out the
% possibility of genotyping errors.

% First merge the two control cohorts.
clear
fprintf('Merging controls.\n');
a  = load('wtccc_58bc.mat');
b  = load('wtccc_ukbs.mat');
na = size(a.X,1);
nb = size(b.X,1);
[X y labels chr pos minor major] = ...
    mergedata(a.X,zeros(na,1),a.labels,a.chr,a.pos,a.minor,a.major,...
	      b.X,zeros(nb,1),b.labels,b.chr,b.pos,b.minor,b.major);
save('wtccc_controls.mat','X','labels','chr','pos','minor','major','-v7.3');

% BIPOLAR DISORDER.
% I do not filter out any SNPs from the bipolar disorder data set.
clear
fprintf('Creating bipolar disorder data set.\n');
cases    = load('wtccc_bd.mat');
ctrls    = load('wtccc_controls.mat');
numcases = size(cases.X,1);
numctrls = size(ctrls.X,1);
[X y labels chr pos minor major] = ...
    mergedata(cases.X,ones(numcases,1),cases.labels,cases.chr,...
	      cases.pos,cases.minor,cases.major,...
	      ctrls.X,zeros(numctrls,1),ctrls.labels,ctrls.chr,...
	      ctrls.pos,ctrls.minor,ctrls.major);
snps = find(sum(X) < 1);
[X labels chr pos minor major] = removesnps(X,labels,chr,pos,minor,major,snps);
save('bd.mat','X','y','labels','chr','pos','minor','major','-v7.3');

% CORONARY ARTERY DISEASE.  
% I filter out one SNP from the coronary artery disease data set: rs6553488
% on chromosome 4 at position 171.38 Mb, with Bayes factor 1448 and minor
% allele frequency 0.46. No nearby SNP has a single-SNP Bayes factor greater
% than 11.
clear
fprintf('Creating coronary artery disease data set.\n');
cases    = load('wtccc_cad.mat');
ctrls    = load('wtccc_controls.mat');
numcases = size(cases.X,1);
numctrls = size(ctrls.X,1);
[X y labels chr pos minor major] = ...
    mergedata(cases.X,ones(numcases,1),cases.labels,cases.chr,...
	      cases.pos,cases.minor,cases.major,...
	      ctrls.X,zeros(numctrls,1),ctrls.labels,ctrls.chr,...
	      ctrls.pos,ctrls.minor,ctrls.major);
badsnp = 6553488;
snps   = [ find(labels == badsnp) find(sum(X) < 1) ];
[X labels chr pos minor major] = removesnps(X,labels,chr,pos,minor,major,snps);
save('cad.mat','X','y','labels','chr','pos','minor','major','-v7.3');

% CROHN'S DISEASE.
% Based on the quality control criterion, I filter out 2 additional SNPS in
% the Crohn's disease data set: rs1914328 on chromosome 8 at position 69.45
% Mb (BF = 6630, MAF = 0.43) and rs6601764 on chromosome 10 at position 3.85
% Mb (BF = 4277, MAF = 0.43). No nearby SNPs have a single-SNP Bayes factor
% greater than 46.
clear
fprintf('Creating Crohn''s disease data set.\n');
cases    = load('wtccc_cd.mat');
ctrls    = load('wtccc_controls.mat');
numcases = size(cases.X,1);
numctrls = size(ctrls.X,1);
[X y labels chr pos minor major] = ...
    mergedata(cases.X,ones(numcases,1),cases.labels,cases.chr,...
	      cases.pos,cases.minor,cases.major,...
	      ctrls.X,zeros(numctrls,1),ctrls.labels,ctrls.chr,...
	      ctrls.pos,ctrls.minor,ctrls.major);
badsnps    = [ 1914328 6601764 ];
[ans snps] = intersect(labels,badsnps);
snps       = [ snps find(sum(X) < 1) ];
[X labels chr pos minor major] = removesnps(X,labels,chr,pos,minor,major,snps);
save('cd.mat','X','y','labels','chr','pos','minor','major','-v7.3');

% HYPERTENSION.
% I do not filter out any SNPs from the hypertension data set.
clear
fprintf('Creating hypertension data set.\n');
cases    = load('wtccc_ht.mat');
ctrls    = load('wtccc_controls.mat');
numcases = size(cases.X,1);
numctrls = size(ctrls.X,1);
[X y labels chr pos minor major] = ...
    mergedata(cases.X,ones(numcases,1),cases.labels,cases.chr,...
	      cases.pos,cases.minor,cases.major,...
	      ctrls.X,zeros(numctrls,1),ctrls.labels,ctrls.chr,...
	      ctrls.pos,ctrls.minor,ctrls.major);
snps = find(sum(X) < 1);
[X labels chr pos minor major] = removesnps(X,labels,chr,pos,minor,major,snps);
save('ht.mat','X','y','labels','chr','pos','minor','major','-v7.3');

% RHEUMATOID ARTHRITIS
% I do not filter out any SNPs from the rheumatoid arthritis data set.
clear
fprintf('Creating rheumatoid arthritis data set.\n');
cases    = load('wtccc_ra.mat');
ctrls    = load('wtccc_controls.mat');
numcases = size(cases.X,1);
numctrls = size(ctrls.X,1);
[X y labels chr pos minor major] = ...
    mergedata(cases.X,ones(numcases,1),cases.labels,cases.chr,...
	      cases.pos,cases.minor,cases.major,...
	      ctrls.X,zeros(numctrls,1),ctrls.labels,ctrls.chr,...
	      ctrls.pos,ctrls.minor,ctrls.major);
snps = find(sum(X) < 1);
[X labels chr pos minor major] = removesnps(X,labels,chr,pos,minor,major,snps);
save('ra.mat','X','y','labels','chr','pos','minor','major','-v7.3');

% TYPE 1 DIABETES.
% I do not filter out any SNPs from the type 1 diabetes data set.
clear
fprintf('Creating type 1 diabetes data set.\n');
cases    = load('wtccc_t1d.mat');
ctrls    = load('wtccc_controls.mat');
numcases = size(cases.X,1);
numctrls = size(ctrls.X,1);
[X y labels chr pos minor major] = ...
    mergedata(cases.X,ones(numcases,1),cases.labels,cases.chr,...
	      cases.pos,cases.minor,cases.major,...
	      ctrls.X,zeros(numctrls,1),ctrls.labels,ctrls.chr,...
	      ctrls.pos,ctrls.minor,ctrls.major);
snps = find(sum(X) < 1);
[X labels chr pos minor major] = removesnps(X,labels,chr,pos,minor,major,snps);
save('t1d.mat','X','y','labels','chr','pos','minor','major','-v7.3');

% TYPE 2 DIABETES.
% I do not filter out any SNPs from the type 2 diabetes data set.
clear
fprintf('Creating type 2 diabetes data set.\n');
cases    = load('wtccc_t2d.mat');
ctrls    = load('wtccc_controls.mat');
numcases = size(cases.X,1);
numctrls = size(ctrls.X,1);
[X y labels chr pos minor major] = ...
    mergedata(cases.X,ones(numcases,1),cases.labels,cases.chr,...
	      cases.pos,cases.minor,cases.major,...
	      ctrls.X,zeros(numctrls,1),ctrls.labels,ctrls.chr,...
	      ctrls.pos,ctrls.minor,ctrls.major);
snps = find(sum(X) < 1);
[X labels chr pos minor major] = removesnps(X,labels,chr,pos,minor,major,snps);
save('t2d.mat','X','y','labels','chr','pos','minor','major','-v7.3');
