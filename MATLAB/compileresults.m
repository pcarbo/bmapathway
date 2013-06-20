% This is a script to compile and display the results from my analysis of
% the case-control studies for all seven diseases.
clear

% Which part of the results to show.
part = 'A';

% Directories where the data and results are stored.
datadir    = '/tmp/pcarbo';
resultsdir = '/tmp/pcarbo/results';

% Distance for assigning SNPs to genes (da), and maximum number of SNP-gene
% annotations.
da = 100e3;  
na = 1e6;    

% The columns of this table give, from left to right: (1) the
% abbreviation for the disease; (2) the full disease name; (3) colour
% uses to denote the disease in scatterplots and barplots;
d = { 'bd'  'bipolar disorder'        'gold'        'bd/pathway-bd-B.mat'
      'cad' 'coronary artery disease' 'salmon'      'cad/pathway-cad-B.mat'
      'cd'  'Crohn''s disease'        'firebrick'   'cd/pathway-cd-B.mat' 
      'ht'  'hypertension'            'forestgreen' 'ht/pathway-ht-B.mat'
      'ra'  'rheumatoid arthritis'    'skyblue'     'ra/pathway-ra-D.mat' 
      't1d' 'type 1 diabetes'         'royalblue'   't1d/pathway-t1d-D.mat'
      't2d' 'type 2 diabetes'         'black'       't2d/pathway-t2d-B.mat'};
disease     = d(:,1);
diseasename = d(:,2);
colours     = d(:,3);
nullfile    = d(:,4);
clear d

% Load the gene and pathway data.
load(strcat(datadir,'/gene_35.1.mat'));
load(strcat(datadir,'/pathway.mat'));

% Assign shorter labels to some of the pathways.
pathway.nickname       = cell(length(pathway.label),1);
pathway.nickname{125}  = 'Striated muscle contraction';
pathway.nickname{126}  = 'Muscle contraction';
pathway.nickname{214}  = 'S6K1 signaling (PC)';
pathway.nickname{215}  = 'Release of eIF4E';
pathway.nickname{323}  = 'Cyclin D events in G1 (PC)';
pathway.nickname{381}  = 'Chemokine receptors';
pathway.nickname{447}  = 'Immune system (PC)';
pathway.nickname{635}  = 'Oligomerization of connexins';
pathway.nickname{637}  = 'Transport of connexons';
pathway.nickname{639}  = 'Gap jct trafficking';
pathway.nickname{640}  = 'Gap jct assembly';
pathway.nickname{684}  = 'Methionine salvage pathway';
pathway.nickname{705}  = 'Incretin synthesis (PC)';
pathway.nickname{706}  = 'Regulation of GLP-1 (PC)';
pathway.nickname{1012} = 'Aspartate biosynthesis';
pathway.nickname{1136} = 'IPn biosynthesis';
pathway.nickname{1162} = 'Aspartate degradation II';
pathway.nickname{1243} = 'Targets of TAp63';
pathway.nickname{1254} = 'Targets of C-MYC (PC)';
pathway.nickname{1273} = 'IL-27 signaling';
pathway.nickname{1312} = 'CN-regulated transcription (PC)';
pathway.nickname{1343} = 'p38-alpha + beta regulation';
pathway.nickname{1354} = 'p38-MAPK signaling';
pathway.nickname{1360} = 'IL-2 signaling (PC)';
pathway.nickname{1390} = 'Targets of Fra1 and Fra2';
pathway.nickname{1421} = 'IL-23 signaling (PC)';
pathway.nickname{1426} = 'IL-12 signaling';
pathway.nickname{1447} = 'Wnt';
pathway.nickname{1508} = 'Cell cycle G1/S';
pathway.nickname{1589} = 'IL-12 and Stat4 in Th1 cells';
pathway.nickname{1611} = 'Malate-aspartate shuttle';
pathway.nickname{1635} = 'NO2-dependent IL12 in NK cells';
pathway.nickname{1700} = 'Chemokine receptors in T-cells';
pathway.nickname{1727} = 'Th1/Th2 differentiation';
pathway.nickname{1748} = 'Arf inhibits ribosomal RNA';
pathway.nickname{1840} = 'Alanine biosynthesis';
pathway.nickname{1877} = 'Asn and Asp biosynthesis';
pathway.nickname{2069} = 'Cytokine signaling';
pathway.nickname{2088} = 'Measles';
pathway.nickname{2118} = 'TP53 network';
pathway.nickname{2141} = 'Id signaling';
pathway.nickname{2208} = 'Wnt signaling';
pathway.nickname{2222} = 'Ala and Asp metabolism';
pathway.nickname{2288} = 'Regulation of GLP-1 (BS)';
pathway.nickname{2289} = 'Incretin synthesis (BS)';
pathway.nickname{2292} = 'IFN-gamma signaling';
pathway.nickname{2294} = 'IFN-alpha/beta signaling';
pathway.nickname{2326} = 'Targets of C-MYC (BS)';
pathway.nickname{2423} = 'IL signaling';
pathway.nickname{2576} = 'ErbB receptor signaling';
pathway.nickname{2591} = 'IL-23 signaling (BS)';
pathway.nickname{2597} = 'CN-regulated transcription (BS)';
pathway.nickname{2608} = 'IL-2 signaling (BS)';
pathway.nickname{2653} = 'IL-12 signaling (BS)';
pathway.nickname{2827} = 'Immune system (BS)';
pathway.nickname{2795} = 'S6K1 signaling (BS)';
pathway.nickname{3084} = 'Cyclin D events in G1 (BS)';
pathway.nickname{3135} = 'Cys and Met metabolism';
pathway.nickname{3146} = 'Allograft rejection';
pathway.nickname{3149} = 'Asthma';
pathway.nickname{3205} = 'TGF-beta signaling';
pathway.nickname{3240} = 'Sulfur metabolism';
pathway.nickname{3287} = 'Phe, Tyr and Trp biosynthesis';  
pathway.nickname{3322} = 'MHC';
pathway.nickname{3323} = 'xMHC';

switch part
 case 'A'

  % POSTERIOR ESTIMATES OF THETA0
  % -----------------------------
  % Get the posterior mean and 95% credible interval for the genome-wide
  % log-odds (theta0) for each disease.
  n = length(disease);
  E = zeros(n,1);
  a = zeros(n,1);
  b = zeros(n,1);

  % Repeat for each disease.
  for i = 1:n
    
    % Load the results from the analysis.
    fprintf('Loading results for %s.\n',diseasename{i});
    load(strcat(resultsdir,'/',nullfile{i}));

    % Get the normalized importance weights.
    w = normalizelogweights(logw);

    % Compute the posterior mean of the genome-wide log-odds.
    E(i) = dot(theta0,w);

    % Compute the 95% credible interval for the genome-wide log-odds.
    [a(i) b(i)] = cred(theta0,w,E(i));
  end

  % Display posterior estimates of the genome-wide log-odds for each disease. 
  figure(1);
  clf
  windowpos = get(gcf,'Position');
  set(gcf,'Name','Results (Part A)','Color','white','NumberTitle','off',...
          'PaperPositionMode','auto','Position',[windowpos(1:2) 550 400]);
  herrorbar(E,1:n,E-a,b-E,'ko');
  set(gca,'FontSize',10,'FontName','fixed');
  set(gca,'TickDir','out','TickLength',[0.02 0.02]);
  set(gca,'XLim',[-6.5,-2.5],'XTick',-6:-3);
  set(gca,'YDir','reverse','YTick',1:n,'YTickLabel',disease);
  set(gca,'XGrid','on');
  xlabel('genome-wide log-odds')
  ylabel('disease');

 case 'B'

  % REGIONS RELEVANT TO DISEASE UNDER NULL
  % --------------------------------------
  % List selected regions of the genome with moderate to strong evidence for
  % disease risk factors. Specifically, show overlapping 50-SNP segments for
  % which there is at least a 0.5 posterior probability that one or more
  % SNPs in the segment are included in the multi-marker disease model (that
  % is, P1 > 0.5).

  % Show the table column headers.
  fprintf('RESULTS (PART B)\n');
  fprintf('                                            ');
  fprintf('                          MAF   MAF \n');
  fprintf('dis chr segment(Mb)     P1   P2 SNP         ');
  fprintf('PIP A/a   LOR(95%% CI)    ctrls cases\n');

  % Repeat for each disease.
  for i = 1:length(disease)

    % Load the genotype and phenotype data.
    load(strcat(datadir,'/',disease{i},'.mat'));

    % Get the cases and controls.
    cases = find(y == 1);
    ctrls = find(y == 0);

    % Create overlapping segments with 50 SNPs.
    [Aseg segchr segpos] = makesegments2(chr,pos,25,5e4);
    [Aseg segchr segpos] = combinesegments(Aseg,segchr,segpos);

    % Load the results from the analysis.
    load(strcat(resultsdir,'/',nullfile{i}));

    % Get the normalized importance weights.
    w = normalizelogweights(logw);

    % Compute the posterior inclusion probabilities (PIPs) averaged over
    % the settings of the hyperparameters.
    PIP = alpha * w;

    % For each segment, compute the statistics P1 and P2 averaged over
    % settings of the hyperparameters.
    [P1 P2] = compilesegstats(Aseg,alpha,w);

    % Show the segments for which the probability of containing an included
    % SNPs exceeds the specified threshold.
    segs = find(P1 > 0.5)';
    for j = segs
          
      % Get the SNPs in the segment.
      snps = find(Aseg(:,j));
          
      % Get the SNP in the segment with the largest posterior inclusion
      % probability.
      [ans k] = max(PIP(snps));
      k       = snps(k);
      
      % Get the mean and 95% credible interval for the additive effect
      % (i.e. log-odds ratio) corresponding to the SNP with the largest
      % PIP in the segment.
      E     = dot(w,mu(k,:));
      beta  = (-2:0.01:2)';
      r     = normpdf(repmat(beta,1,length(w)),...
                      repmat(mu(k,:),length(beta),1),...
                      repmat(s(k,:),length(beta),1));
      [a b] = cred(beta,r*w,E);
      
      % Print statistics about the selected segment.
      str = sprintf('%0.2f-%-0.2f',min(pos(snps))/1e6,max(pos(snps))/1e6);
      fprintf('%3s %3d %-13s %0.2f %0.2f %-10s %0.2f %c/%c ',...
              upper(disease{i}),segchr(j),str,P1(j),P2(j),...
              sprintf('rs%d',labels(k)),PIP(k),major(k),minor(k));
      fprintf('%+0.2f(%0.2f-%0.2f) %0.3f %0.3f\n',...
              E,min(abs([a b])),max(abs([a b])),...
              maf(X(ctrls,k)),maf(X(cases,k)));
    end

    % Delete the genotype matrix so that we have enough space in memory
    % to load the data for the next disease.
    clear X
  end

 case 'C'

  % GENOME-WIDE SCANS UNDER NULL
  % ----------------------------
  % Show the genome-wide scan (the "Manhattan plot") for all seven
  % diseases. 
  figure(1);
  clf
  windowpos = get(gcf,'Position');
  set(gcf,'Name','Results (Part C)','Color','white','NumberTitle','off',...
          'PaperPositionMode','auto','Position',[windowpos(1:2) 675 850]);

  % Repeat for each disease.
  n = length(disease);
  for i = 1:n

    % Load the genotype and phenotype data.
    fprintf('Loading data for %s.\n',diseasename{i});
    load(strcat(datadir,'/',disease{i},'.mat'));

    % Create overlapping segments with 50 SNPs.
    fprintf('Partitioning the genome into segments.\n');
    [Aseg segchr segpos] = makesegments2(chr,pos,25,5e4);
    [Aseg segchr segpos] = combinesegments(Aseg,segchr,segpos);

    % Load the results from the analysis.
    fprintf('Loading results for %s.\n',diseasename{i});
    load(strcat(resultsdir,'/',nullfile{i}));

    % Get the normalized importance weights.
    w = normalizelogweights(logw);

    % For each segment, compute the statistics P1 and P2 averaged over
    % settings of the hyperparameters.
    [P1 P2] = compilesegstats(Aseg,alpha,w);

    % Show the genome-wide scan, highlighting the regions of the genome
    % that are listed in the table (Part B) because they have moderate to
    % strong evidence for disease risk factors.
    segs = find(P1 > 0.5);
    subplot(n,1,i);
    genomewideplot(segchr,segpos,P1,segs);
    set(gca,'FontSize',9,'FontName','fixed');
    set(gca,'YLim',[0 1],'YTick',[0 0.5 1]);
    set(gca,'TickLength',[0.01 0.01]);
    ylabel('PIP');
    title(diseasename{i});

    % Delete the genotype matrix so that we have enough space in memory
    % to load the data for the next disease.
    clear X
  end

 case 'D'

  % P1 SCATTERPLOT FOR RA AND T1D
  % -----------------------------
  % For rheumatoid arthritis and type 1 diabetes, draw a scatterplot comparing
  % the posterior statistics (P1) from Stage C, where we analyze all SNPs
  % across the genome simultaneously, and Stage D, where we treat SNPs in
  % the the MHC region separately from SNPs outside the MHC.
  figure(1);
  clf
  windowpos = get(gcf,'Position');
  set(gcf,'Name','Results (Part D)','Color','white','NumberTitle','off',...
          'PaperPositionMode','auto','Position',[windowpos(1:2) 975 375]);

  % Repeat for each disease under consideration (here I'm just looking at
  % rheumatoid arthritis and type 1 diabetes).
  j = 0;
  for i = [5 6]
    j = j + 1;

    % Load the genotype and phenotype data.
    fprintf('Loading data for %s.\n',diseasename{i});
    load(strcat(datadir,'/',disease{i},'.mat'));

    % Create overlapping segments with 50 SNPs.
    fprintf('Partitioning the genome into segments.\n');
    [Aseg segchr segpos] = makesegments2(chr,pos,25,5e4);
    [Aseg segchr segpos] = combinesegments(Aseg,segchr,segpos);

    % Determine which of the segments lie in the ("extended") MHC region.
    isinmhc = (Aseg' * inmhcregion(chr,pos)) > 0;
    
    % Load the results from Stages C and D of the analysis.
    fprintf('Loading results for %s.\n',diseasename{i});
    c = load(sprintf('%s/%s/pathway-%s-C.mat',...
                     resultsdir,disease{i},disease{i}));
    d = load(sprintf('%s/%s/pathway-%s-D.mat',...
                     resultsdir,disease{i},disease{i}));

    % Get the normalized importance weights.
    c.w = normalizelogweights(c.logw);
    d.w = normalizelogweights(d.logw);

    % For each segment, compute the statistics P1 and P2 averaged over
    % settings of the hyperparameters.
    [c.P1 c.P2] = compilesegstats(Aseg,c.alpha,c.w);
    [d.P1 d.P2] = compilesegstats(Aseg,d.alpha,d.w);

    % Show scatterplot of P1 from Stages C and D.
    subplot(1,2,j);
    plot(c.P1(isinmhc),d.P1(isinmhc),'.',...
         'Color',rgb('firebrick'),'MarkerSize',20);
    hold on
    plot(c.P1(~isinmhc),d.P1(~isinmhc),'.',...
         'Color',rgb('dodgerblue'),'MarkerSize',20);
    hold off
    box on
    set(gca,'FontSize',10,'FontName','fixed');
    set(gca,'XLim',[-0.05 1.05],'XTick',[0 0.5 1]);
    set(gca,'YLim',[-0.05 1.05],'YTick',[0 0.5 1]);
    set(gca,'TickDir','out','TickLength',[0.015 0.015]);
    xlabel('P1 when all SNPs are analyzed jointly');
    ylabel('P1 when MHC is analyzed separately');
    title(diseasename{i});

    % Label some of the segments by the chromosome and position on the
    % chromosome.
    I = find(c.P1 > 0.25 | d.P1 > 0.25)';
    for t = I
      text(c.P1(t),d.P1(t),sprintf(' %d %0.2f',segchr(t),segpos(t)/1e6),...
           'VerticalAlignment','bottom','FontSize',9,'FontName','fixed');
    end

    % Delete the genotype matrix so that we have enough space in memory
    % to load the data for the next disease.
    clear X
  end

 case 'E'

  % SINGLE-MARKER VS MULTI-MARKER SCATTERPLOTS
  % ------------------------------------------
  % Load the table of genome-wide associations for all 7 diseases.
  a = importdata('results/fourthdraft/p1vspvalue.txt');
  t = a.textdata(2:end,:);
  a = a.data;

  % Get the columns of the table.
  dis     = t(:,1);
  locus   = t(:,2);
  chr     = t(:,3);
  region  = t(:,4);
  P1      = a(:,1);
  trend   = a(:,2);
  logBF   = a(:,3);

  % Cut off the trend p-values so they all lie between 1e-5 and 1e-10, and
  % cut off the (log) Bayes factors so that they are at most 9 and at
  % least 3.
  trend = min(1e-5,max(1e-10,trend));
  logBF = max(3,min(9,logBF));

  % Jitter the top-right points in the scatterplot.
  I        = find(P1 == 1 & logBF == 9);
  n        = length(I);
  P1(I)    = P1(I) + 0.02*randn(n,1);
  trend(I) = 10.^(log10(trend(I)) + 0.1*randn(n,1));
  logBF(I) = logBF(I) + 0.1*randn(n,1);

  figure(1);
  clf
  windowpos = get(gcf,'Position');
  set(gcf,'Name','Results (Part E)','Color','white','NumberTitle','off',...
          'PaperPositionMode','auto','Position',[windowpos(1:2) 1200 450]);
  
  % Show scatterplot of trend p-value vs P1.
  subplot(1,2,1)
  hold on
  for i = 1:length(disease)
    I = find(strcmp(upper(disease{i}),dis));
    plot(log10(trend(I)),P1(I),'.','Color',rgb(colours{i}),'MarkerSize',20);
  end

  % Highlight some of the regions in the scatterplot.
  I = find(~strcmp(locus,'-'))';
  for i = I
    text(log10(trend(i)),P1(i),sprintf(' %s',locus{i}),...
         'VerticalAlignment','bottom','FontSize',10,'FontName','fixed');
  end
  hold off
  box on
  set(gca,'FontSize',10,'FontName','fixed');
  set(gca,'XLim',[-10.5 -4.5],'XTick',[-10 -8 -6.3 -5],'XDir','reverse');
  set(gca,'YLim',[-0.1 1.1],'YTick',[0 0.5 1]);
  set(gca,'TickDir','out','TickLength',[0.015 0.015]);
  xlabel('log_{10} trend p-value');
  ylabel('P1');

  % Show scatterplot of trend p-value vs (additive) Bayes factor.
  subplot(1,2,2)
  hold on
  for i = 1:length(disease)
    I = find(strcmp(upper(disease{i}),dis) & logBF >= 0);
    plot(logBF(I),P1(I),'.','Color',rgb(colours{i}),'MarkerSize',20);
  end

  % Highlight some of the regions in the scatterplot.
  I = find(~strcmp(locus,'-'))';
  for i = I
    text(logBF(i),P1(i),sprintf(' %s',locus{i}),...
         'VerticalAlignment','bottom','FontSize',10,'FontName','fixed');
  end
  hold off
  box on
  legend(disease,'Location','SouthEast');
  legend boxoff
  set(gca,'FontSize',10,'FontName','fixed');
  set(gca,'XLim',[2.5 9.5],'XTick',3:3:9);
  set(gca,'YLim',[-0.1 1.1],'YTick',[0 0.5 1]);
  set(gca,'TickDir','out','TickLength',[0.015 0.015]);
  xlabel('log_{10} Bayes factor');
  ylabel('P1');

 case 'F'

  % PATHWAY STATISTICS
  % ------------------ 
  % Get the set of pathways included in the analysis. We ignore pathways with
  % less than one gene assigned to the reference genome. We also ignore
  % NCI/Nature PID pathways from Pathway Commons with more than 500 genes.
  numgenes = pathway.genes' * inref(gene);
  bad      = strcmp(pathway.database,'PC') & ...
             strcmp(pathway.source,'PID') & ...
             numgenes > 500;
  analyzed = numgenes > 1 & ~bad;
  if any(analyzed ~= full(sumcols(getinitpaths(gene,pathway))))
    error('We miscalculated which pathways are included in the analysis')
  end
  paths = find(numgenes < 2);
  fprintf('%d pathways have less than 2 genes mapped to ',...
          length(paths) + sum(cellfun(@numel,pathway.synonyms(paths))));
  fprintf('reference genome.\n');
  paths = find(bad);
  fprintf('%d pathways from NCI Nature PID (Pathway Commons) ',...
          length(paths) + sum(cellfun(@numel,pathway.synonyms(paths))));
  fprintf('are omitted.\n\n');

  % Show the number of pathways (unique gene sets) from each database
  % that are included in the analysis.
  sources   = { 'Reactome'; 'KEGG'; 'BioCarta'; 'HumanCyc'; 'PID';
                'WikiPathways'; 'Panther'; 'Cell Map'; 'NA' };
  databases = { 'PC'; 'BioSystems'; 'BioCarta'; 'Panther'; 'NA' };

  % Show the table column headers.
  fprintf('RESULTS (PART D)\n');
  fprintf('database     source     #sets genes\n');

  % Repeat for each source (e.g. Reactome), and for each database
  % (e.g. Pathway Commons).
  covered  = zeros(length(gene.symbol),1);
  numpaths = zeros(length(databases),length(sources));
  numgenes = zeros(length(databases),length(sources));
  for i = 1:length(sources)
    for j = 1:length(databases)

      % Get the number of pathways (unique gene sets) from the database
      % and source that are included in the analysis.
      paths = find(strcmp(sources{i},pathway.source) & ...
                   strcmp(databases{j},pathway.database) & ...
                   analyzed);
      numpaths(j,i) = length(paths);

      % Get the number of additional genes that are assigned to at least
      % one of these pathways.
      genes = find(~covered & inref(gene) & ...
                   sumcols(pathway.genes(:,paths)));
      numgenes(j,i) = length(genes);
      covered(genes) = 1;

      if numpaths(j,i) > 0
        fprintf('%-12s %-10s %5d %5d\n',sources{i},databases{j},...
                numpaths(j,i),numgenes(j,i));
      end
    end
  end
  totalnumgenes = sum(inref(gene));
  fprintf('total                   %5d %5d out of %d\n',...
          sum(analyzed),sum(covered),totalnumgenes);

  % Show the same results in a different form using a bar chart.
  numpaths = numpaths(numpaths > 100);
  numpaths = [ numpaths
               sum(analyzed) - sum(numpaths) ];
  numgenes = numgenes(numgenes > 1e3);
  numgenes = [ numgenes
               sum(covered) - sum(numgenes) ];

  % Show the same results in a different form using a bar chart.
  figure(1);
  clf
  windowpos = get(gcf,'Position');
  set(gcf,'Name','Results (Part F)','Color','white','NumberTitle','off',...
          'PaperPositionMode','auto','Position',[windowpos(1:2) 550 400]);

  subplot(1,2,1);
  bar([numpaths zeros(size(numpaths))]','BarLayout','stacked',...
      'BarWidth',0.5,'EdgeColor',rgb('white'),'FaceColor',rgb('silver'),...
      'LineWidth',2);
  set(gca,'FontSize',10,'FontName','fixed');
  set(gca,'XLim',[0 2],'XTick',[]);
  set(gca,'YLim',[0 sum(analyzed)],'YTick',[]);
  ylabel('number of pathways');  
  t = cumsum(numpaths);
  text(repmat(0.7,size(t)),t,num2str(numpaths),'HorizontalAlignment',...
       'right','VerticalAlignment','Top','FontSize',9,'FontName','fixed');

  subplot(1,2,2);
  bar(totalnumgenes,'BarWidth',0.5,'EdgeColor',rgb('white'),...
      'FaceColor',rgb('orange'),'LineWidth',2);
  hold on
  bar([numgenes zeros(size(numgenes))]','BarLayout','stacked',...
      'BarWidth',0.5,'EdgeColor',rgb('white'),'FaceColor',rgb('silver'),...
      'LineWidth',2);
  hold off
  set(gca,'FontSize',10,'FontName','fixed');
  set(gca,'XLim',[0 2],'XTick',[]);
  set(gca,'YLim',[0 totalnumgenes],'YTick',[0 sum(covered) totalnumgenes],...
          'YTickLabel',[0 sum(covered) totalnumgenes]);
  ylabel('gene coverage');
  t = cumsum(numgenes);
  text(repmat(0.7,size(t)),t,num2str(numgenes),'HorizontalAlignment',...
       'right','VerticalAlignment','Top','FontSize',9,'FontName','fixed');

 case 'G'

  % MORE PATHWAY STATISTICS
  % -----------------------
  % Load the Crohn's disease data.
  fprintf('Loading Crohn''s disease data.\n');
  load(strcat(datadir,'/cd.mat'));
  
  % ASSIGN SNPS TO GENES, THEN GENES TO PATHWAYS.
  fprintf('Assigning SNPs to genes and pathways.\n');
  Agene  = snps2genes(chr,pos,gene,da,na);
  A      = spones(Agene * pathway.genes);
  A      = modifymhc(pathway,A);

  % Get the set of pathways included in the analysis.
  numgenes = pathway.genes' * inref(gene);
  analyzed = full(sumcols(getinitpaths(gene,pathway)) > 0);

  % Show statistics about assignment of genes and SNPs to pathways.
  figure(1);
  clf
  windowpos = get(gcf,'Position');
  set(gcf,'Name','Results (Part G)','Color','white','NumberTitle','off',...
          'PaperPositionMode','auto','Position',[windowpos(1:2) 550 425]);

  % Show the histogram of gene set sizes.
  subplot(2,2,1);
  edges  = [ 2 10 100 1e3 1e4 ];
  n      = length(edges);
  counts = histc(numgenes(analyzed),edges);
  counts = counts(1:n-1);
  h      = barh(counts,0.7,'histc');
  text(counts + 50,0.5 + (1:n-1),num2str(counts,'%-d'),...
       'FontSize',10,'FontName','fixed');
  set(h,'LineStyle','none','FaceColor',rgb('silver'));
  set(gca,'FontSize',10,'FontName','fixed');
  set(gca,'XLim',[0 2400],'XTick',[]);
  set(gca,'YLim',[0 n+1],'YDir','reverse','YTick',1:n,'YTickLabel',edges);
  set(gca,'TickLength',[0 0]);
  xlabel('pathways');
  ylabel('genes');

  % Show the histogram of number of SNPs assigned to pathways.
  subplot(2,2,2);
  edges   = [ 0 10 100 1e3 1e4 1e5 ];
  n       = length(edges);
  numsnps = full(sumrows(A))';
  counts  = histc(numsnps(analyzed),edges);
  counts  = counts(1:n-1);
  h       = barh(counts,0.7,'histc');
  text(counts + 50,0.5 + (1:n-1),num2str(counts,'%-d'),...
       'FontSize',10,'FontName','fixed');
  set(h,'LineStyle','none','FaceColor',rgb('silver'));
  set(gca,'FontSize',10,'FontName','fixed');
  set(gca,'XLim',[0 2200],'XTick',[]);
  set(gca,'YLim',0.5 + [0 n],'YDir','reverse','YTick',1:n,'YTickLabel',edges);
  set(gca,'TickLength',[0 0]);
  xlabel('pathways');
  ylabel('SNPs');

  % Show the histogram of number of pathways assigned to each gene.
  subplot(2,2,3);
  edges    = [ 0 1 10 100 1e3 ];
  n        = length(edges);
  numpaths = full(sumcols(pathway.genes(inref(gene),analyzed)));
  counts   = histc(numpaths,edges);
  counts   = counts(1:n-1);
  h        = barh(counts,0.7,'histc');
  text(counts + 500,0.5 + (1:n-1),num2str(counts,'%-d'),...
       'FontSize',10,'FontName','fixed');
  set(h,'LineStyle','none','FaceColor',rgb('silver'));
  set(gca,'FontSize',10,'FontName','fixed');
  set(gca,'XLim',[0 2e4],'XTick',[]);
  set(gca,'YLim',[0 n+1],'YDir','reverse','YTick',1:n,'YTickLabel',edges);
  set(gca,'TickLength',[0 0]);
  xlabel('genes');
  ylabel('pathways');

  % Show the histogram of number of pathways assigned to each SNP.
  subplot(2,2,4);
  edges    = [ 0 1 10 100 1e3 ];
  n        = length(edges);
  numpaths = full(sumcols(A(:,analyzed)));
  counts   = histc(numpaths,edges);
  counts   = counts(1:n-1);
  h        = barh(counts,0.7,'histc');
  text(counts + 1e4,0.5 + (1:n-1),num2str(counts,'%-d'),...
     'FontSize',10,'FontName','fixed');
  set(h,'LineStyle','none','FaceColor',rgb('silver'));
  set(gca,'FontSize',10,'FontName','fixed');
  set(gca,'XLim',[0 3.5e5],'XTick',[]);
  set(gca,'YLim',[0 n+1],'YDir','reverse','YTick',1:n,'YTickLabel',edges);
  set(gca,'TickLength',[0 0]);
  xlabel('SNPs');
  ylabel('pathways');

 case 'H'

  % DISTRIBUTION OF BAYES FACTORS
  % -----------------------------
  % The columns of this table specify (1) the disease; (2) the colour used to
  % draw the bars in the bar plots; (3) the file storing the Bayes factor
  % estimates.
  D = { 'bd'  'dodgerblue' 'bd/pathway-bd-C.mat'
        'cad' 'dodgerblue' 'cad/pathway-cad-C.mat'
        'cd'  'dodgerblue' 'cd/pathway-cd-C.mat'
        'ht'  'dodgerblue' 'ht/pathway-ht-C.mat'
        'ra'  'dodgerblue' 'ra/pathway-ra-E.mat'
        't1d' 'dodgerblue' 't1d/pathway-t1d-E.mat'
        't2d' 'dodgerblue' 't2d/pathway-t2d-C.mat'
        'ra'  'darkorange' 'ra/pathway-ra-I.mat'
        't1d' 'darkorange' 't1d/pathway-t1d-I.mat' };

  d           = D(:,1);
  colours     = D(:,2);
  resultsfile = D(:,3);
  clear D

  figure(1);
  clf
  windowpos = get(gcf,'Position');
  set(gcf,'Name','Results (Part H)','Color','white','NumberTitle','off',...
          'PaperPositionMode','auto','Position',[windowpos(1:2) 675 600]);

  % Repeat for each set of results.
  n = length(d);
  for i = 1:n

    % Find the disease index.
    j = find(strcmp(disease,d{i}));

    % Load the Bayes factors. Set all Bayes factors less than 0.01 to be equal
    % to 0.01, since we are plotting the Bayes factors on the logarithmic
    % scale.
    fprintf('Loading results for %s.\n',diseasename{j});
    load(strcat(resultsdir,'/',resultsfile{i}));
    BF = max(0.01,BF);

    % Show the histogram of the Bayes factors on the logarithmic scale. I
    % label each bar with a number, which is useful for the counts that
    % are very small.
    subplot(3,3,i);
    edges  = [-100 0:4 100];
    ne     = length(edges);
    counts = histc(log10(BF),edges);
    counts = counts(1:ne-1);
    h      = bar(counts,0.7,'histc');
    set(h,'LineStyle','none','FaceColor',rgb(colours{i}));
    set(gca,'FontSize',10,'FontName','fixed');
    set(gca,'XLim',[0.5 ne+0.5],'XTick',1:ne,'XTickLabel',edges);
    set(gca,'YLim',[0 3500],'YTick',[]);
    set(gca,'TickLength',[0 0]);
    xlabel('log_{10}Bayes factor','FontSize',10);
    ylabel('no. gene sets','FontSize',10);
    title(diseasename{j},'FontSize',10);

    % Add the count to the top of each bar in the bar plot.
    for t = 1:ne-1
      if counts(t) > 0
        text(0.5+t,counts(t) + 50,num2str(counts(t)),...
             'FontSize',9,'FontName','fixed',...
             'HorizontalAlignment','center',...
             'VerticalAlignment','bottom');
      end
    end
  end

 case 'I'

  % BAYES FACTOR SCATTERPLOT WITH COARSE AND FINE GRIDS
  % ---------------------------------------------------
  % The columns of this table specify (1) the disease; (2) the file containing
  % rough numerical estimates of the Bayes factors using a coarse grid over
  % the hyperparameters; (3) the file containing more accurate numerical
  % estimates of the Bayes factors using a more fine-grained grid over the
  % hyperparameters.
  %
  % Note that I've omitted the results for RA and T1D that do not condition
  % on enrichment of the MHC and xMHC. Also note that the Bayes factors that
  % *do* condition on enrichment of the MHC are computed with respect to the
  % hypothesis that only the MHC is enriched.
  D = { 'bd'  'bd/pathway-bd-C.mat'   'bd/pathway-bd-D.mat'
        'cad' 'cad/pathway-cad-C.mat' 'cad/pathway-cad-E.mat'
        'cd'  'cd/pathway-cd-C.mat'   'cd/pathway-cd-D.mat'
        'cd'  'cd/pathway-cd-E.mat'   'cd/pathway-cd-G.mat'
        'cd'  'cd/pathway-cd-H.mat'   'cd/pathway-cd-J.mat'
        'ht'  'ht/pathway-ht-C.mat'   'ht/pathway-ht-D.mat' 
        'ra'  'ra/pathway-ra-I.mat'   'ra/pathway-ra-J.mat' 
        'ra'  'ra/pathway-ra-K.mat'   'ra/pathway-ra-L.mat' 
        'ra'  'ra/pathway-ra-M.mat'   'ra/pathway-ra-N.mat' 
        't1d' 't1d/pathway-t1d-I.mat' 't1d/pathway-t1d-J.mat' 
        't1d' 't1d/pathway-t1d-K.mat' 't1d/pathway-t1d-L.mat'  
        't1d' 't1d/pathway-t1d-M.mat' 't1d/pathway-t1d-N.mat' 
        't2d' 't2d/pathway-t2d-C.mat' 't2d/pathway-t2d-E.mat' };

  d            = D(:,1);
  resultsfile1 = D(:,2);
  resultsfile2 = D(:,3);

  figure(1);
  clf
  windowpos = get(gcf,'Position');
  set(gcf,'Name','Results (Part I)','Color','white','NumberTitle','off',...
          'PaperPositionMode','auto','Position',[windowpos(1:2) 550 450]);

  % Plot the x = y line in the scatterplot.
  plot([0 20],[0 20],':','Color',rgb('black'),'LineWidth',1);

  % Repeat for each set of results.
  n = length(d);
  hold on
  for i = 1:n

    % Find the the disease index.
    j = find(strcmp(disease,d{i}));

    % Load the enrichment hypotheses and Bayes factors.
    fprintf('Loading results for %s.\n',diseasename{j});
    a   = load(strcat(resultsdir,'/',resultsfile1{i}));
    H1  = a.H;
    BF1 = a.BF;

    a   = load(strcat(resultsdir,'/',resultsfile2{i}));
    H2  = a.H;
    BF2 = a.BF;
    clear a

    % Find the enrichment hypotheses in the initial computation of Bayes
    % factors that match the enrichment hypotheses in the second attempt
    % at calculating the Bayes factors.
    [ans I] = ismember(H2',H1','rows');

    % Plot the Bayes factors on the logarithmic scale.
    plot(log10(BF1(I)),log10(BF2),'.','Color',rgb(colours{j}),...
         'MarkerSize',20);
  end
  hold off
  box on
  legend(['NA'; d],'Location','SouthEast');
  legend boxoff
  set(gca,'FontSize',10,'FontName','fixed');
  set(gca,'XLim',[0 18],'XTick',0:5:20);
  set(gca,'YLim',[0 18],'YTick',0:5:20);
  set(gca,'TickDir','out','TickLength',[0.01 0.01]);
  xlabel('Rough estimate of log_{10}BF');
  ylabel('More accurate estimate of log_{10}BF');

 case 'J'

  % TABLE OF TOP-RANKED ENRICHED PATHWAYS
  % -------------------------------------
  % The columns of this table specify: (1) the disease; (2) the file
  % containing the more accurate estimates of the Bayes factors; (3) the
  % index of the selected enrichment hypothesis; (4) if equal to 1, the
  % Bayes factor is calculated with respect ot the hypothesis that the MHC
  % is enriched for disease risk factors, and the enrichment hypothesis
  % includes enrichment of the MHC region at a different enrichment rate;
  % (5) if this number is greater than 1, then the Bayes factor is
  % calculated with respect to the hypothesis that the MHC is enriched for
  % disease risk factors, in which case we need to multiply by this number
  % to obtain the Bayes factor with respect to the hypothesis that no
  % pathways are enriched.
  D = { 't1d' 't1d/pathway-t1d-J.mat' 16 1 2.70e54 % IL-2 signaling + MHC.
        'ra'  'ra/pathway-ra-J.mat'   3  1 2.43e21 % Measles + MHC. 
        'cd'  'cd/pathway-cd-D.mat'   11 0 1       % Cytokine signaling.
        't2d' 't2d/pathway-t2d-E.mat' 1  0 1       % Incretin synthesis.
        'cad' 'cad/pathway-cad-E.mat' 4  0 1       % Arf inhibits RNA.
        'bd'  'bd/pathway-bd-D.mat'   2  0 1       % Transport of connetions.
        'ht'  'ht/pathway-ht-D.mat'   11 0 1 };    % Ala biosynthesis.

  disease       = D(:,1);
  resultsfile   = D(:,2);
  hind          = cell2mat(D(:,3));
  mhcisenriched = cell2mat(D(:,4));
  BFmhc         = cell2mat(D(:,5));

  % This is the pathway index corresponding to the major histocompatibility
  % complex (MHC).
  mhcind = 3322;

  figure(1);
  clf
  windowpos = get(gcf,'Position');
  set(gcf,'Name','Results (Part J)','Color','white',...
          'NumberTitle','off','PaperPositionMode','auto',...
          'Position',[windowpos(1:2) 1000 250]);

  % Initialize storage for the results to be shown in the table.
  n           = length(disease);
  paths       = zeros(n,1);
  numgenes    = zeros(n,1);
  numsnps     = zeros(n,1);
  numgenesmhc = zeros(n,1);
  numsnpsmhc  = zeros(n,1);
  E0          = zeros(n,1);
  a0          = zeros(n,1);
  b0          = zeros(n,1);
  E           = zeros(n,1);
  a           = zeros(n,1);
  b           = zeros(n,1);

  % Repeat for each disease. In this loop, I use the following notation:
  % 
  %   i  disease
  %   j  pathway
  %   k  selected enrichment hypothesis
  % 
  for i = 1:n
      
    % Load the Bayes factors and other statistics.
    load(strcat(resultsdir,'/',resultsfile{i}));

    % Select data about the requested enrichment hypothesis.
    k     = hind(i);
    BF    = BF(k);
    H     = H(:,k);
    logw1 = logw1(:,:,k);

    % Get the index of the enriched pathway.
    j        = find(H);
    paths(i) = j;

    % Get the number of genes assigned to the enriched pathway, and the number
    % of genes assigned to the MHC if the MHC region is enriched for disease
    % associations. If the enrichment hypothesis includes enrichment of the
    % MHC, then remove all genes within the MHC from this count.
    if mhcisenriched(i)
      numgenesmhc(i) = sum(full(pathway.genes(:,mhcind) & inref(gene)));
      numgenes(i)    = sum(pathway.genes(:,j) & ~pathway.genes(:,mhcind) & ...
                           inref(gene));
    else
      numgenesmhc(i) = 0;
      numgenes(i)    = sum(pathway.genes(:,j) & inref(gene));
    end

    % Get the number of SNPs assigned to the enriched pathway, and the
    % number of SNPs assigned to the MHC if the MHC region is enriched
    % for disease associations. If the enriched hypothesis includes
    % enrichment of the MHC, then remove all SNPs within the MHC from
    % this count.
    if mhcisenriched(i)
      numsnpsmhc(i) = sum(A(:,mhcind));
      numsnps(i)    = sum(A(:,j) & ~A(:,mhcind));
    else
      numsnpsmhc(i) = 0;
      numsnps(i)    = sum(A(:,j));
    end

    % Get the posterior mean and 95% credible interval interval for the
    % genome-wide log-odds (theta0).
    w             = normalizelogweights(logw1);
    E0(i)         = dot(theta0,sumcols(w));
    [a0(i) b0(i)] = cred(theta0,sumcols(w),E0(i));

    % Get the posterior mean and 95% credible interval interval for the
    % enrichment parameter (theta)
    E(i)        = dot(theta,sumrows(w));
    [a(i) b(i)] = cred(theta,sumrows(w),E(i));

    % Print more details about the enrichment hypothesis.
    if isempty(pathway.synonyms{j})
      star = '';
    else
      star = '*';
    end
    if mhcisenriched(i)
      fprintf('%s: %s + MHC:\n',upper(disease{i}),pathway.nickname{j});
    else
      fprintf('%s: %s\n',upper(disease{i}),pathway.nickname{j});
    end
    fprintf('%s%s (%s/%s,i=%d)\n',pathway.label{j},star,...
            pathway.source{j},pathway.database{j},j);
    BF = BF * BFmhc(i);
    if BF > 1e4
      fprintf('BF     = %0.2e\n',BF);
    else
      fprintf('BF     = %d\n',round(BF));
    end
    fprintf('theta0 = %0.2f (%0.1f,%0.1f)\n',a0(i),E0(i),b0(i));
    fprintf('theta  = %0.2f (%0.1f,%0.1f)\n',a(i),E(i),b(i));
    fprintf('\n');
  end

  % Plot the number of genes.
  subplot(1,6,2:3);
  h = barh(1:n,numgenesmhc + numgenes,0.7);
  set(h,'EdgeColor',rgb('white'),'LineWidth',3,'FaceColor',rgb('darkorange'));
  hold on
  h = barh(1:n,numgenesmhc,0.7);
  set(h,'EdgeColor',rgb('white'),'LineWidth',3,'FaceColor',rgb('darkorange'));
  hold off
  axespos = get(gca,'Position');
  set(gca,'Position',[axespos(1) 0.15 axespos(3) 0.75]);
  set(gca,'FontSize',9,'FontName','fixed');
  set(gca,'YDir','reverse','YLim',[0 n+1],'YTick',1:n);
  set(gca,'XLim',[0 320],'XTick',0:100:300);
  set(gca,'TickLength',[0 0]);
  xlabel('no. genes','FontSize',10);
  str = num2str(numgenes,'%-d');
  I   = find(~mhcisenriched);
  text(numgenes(I) + 10,I,str(I,:),'FontSize',9,'FontName','fixed');
  str = strcat(num2str(numgenesmhc),'+',str);
  I   = find(mhcisenriched);
  text(numgenesmhc(I) + numgenes(I) + 10,I,str(I,:),...
       'FontSize',9,'FontName','fixed');

  % Label the bar plot with the disease and the enriched pathway(s).
  str    = strcat(upper(disease),'-',pathway.nickname(paths));
  I      = find(mhcisenriched);
  str(I) = strcat(str(I),'+','MHC');
  set(gca,'YTickLabel',str);
  
  % Plot the number of SNPs.
  subplot(1,6,4);
  h = barh(1:n,numsnpsmhc + numsnps,0.7);
  set(h,'EdgeColor',rgb('white'),'LineWidth',3,'FaceColor',rgb('darkorange'));
  hold on
  h = barh(1:n,numsnpsmhc,0.7);
  set(h,'EdgeColor',rgb('white'),'LineWidth',3,'FaceColor',rgb('darkorange'));
  hold off
  axespos = get(gca,'Position');
  set(gca,'Position',[axespos(1) 0.15 axespos(3) 0.75]);
  set(gca,'FontSize',9,'FontName','fixed');
  set(gca,'YDir','reverse','YLim',[0 n+1],'YTick',[]);
  set(gca,'XLim',[0 1e4],'XTick',[0 1e4]);
  xlabel('no. SNPs','FontSize',10);
  str = num2str(numsnps,'%-d');
  I   = find(~mhcisenriched);
  text(numsnps(I) + 200,I,str(I,:),'FontSize',9,'FontName','fixed');
  str = strcat(num2str(numsnpsmhc),'+',str);
  I   = find(mhcisenriched);
  text(numsnpsmhc(I) + numsnps(I) + 200,I,str(I,:),...
       'FontSize',9,'FontName','fixed');

  % Show the posterior mean and 95% credible interval of the genome-wide
  % log-odds (theta0) for each enrichment hypothesis.
  subplot(1,6,5);
  h = herrorbar(E0,1:n,E0-a0,b0-E0,'ko');
  set(h,'Markersize',4,'MarkerFaceColor',rgb('black'));
  axespos = get(gca,'Position');
  set(gca,'Position',[axespos(1) 0.15 axespos(3) 0.75]);
  set(gca,'FontSize',9,'FontName','fixed');
  set(gca,'YDir','reverse','YLim',[0 n+1],'YTick',[]);
  set(gca,'XLim',[-6.25 -2],'XTick',-6:-2,'XGrid','on');
  set(gca,'TickLength',[0 0]);
  xlabel('genome-wide log-odds','FontSize',10);

  % Show the posterior mean and 95% credible interval of the enrichment
  % parameter (theta) for each enrichment hypothesis.
  subplot(1,6,6);
  h = herrorbar(E,1:n,E-a,b-E,'ko');
  set(h,'Markersize',4,'MarkerFaceColor',rgb('black'));
  axespos = get(gca,'Position');
  set(gca,'Position',[axespos(1) 0.15 axespos(3) 0.75]);
  set(gca,'FontSize',9,'FontName','fixed');
  set(gca,'YDir','reverse','YLim',[0 n+1],'YTick',[]);
  set(gca,'XLim',[0 5],'XTick',0:5,'XGrid','on');
  set(gca,'TickLength',[0 0]);
  xlabel('log_{10}enrichment','FontSize',10);

 case 'K'

  % EXPANDED RANKING OF ENRICHED PATHWAYS
  % -------------------------------------
  % The columns of this table specify: (1) the disease; (2) the file containing
  % the more accurate estimates of the Bayes factors; (3) if equal to 1, the
  % Bayes factor is calculated with respect ot the hypothesis that the MHC
  % is enriched for disease risk factors, and the enrichment hypothesis
  % includes enrichment of the MHC region at a different enrichment rate;
  % (4) if this number is greater than 1, then the Bayes factor is
  % calculated with respect to the hypothesis that the MHC is enriched for
  % disease risk factors, in which case we need to multiply by this number
  % to obtain the Bayes factor with respect to the hypothesis that no
  % pathways are enriched; (5) the indices of the selected enrichment
  % hypotheses.
  D = { % Top 3 pathways in T1D without conditioning on enrichment of the
        % MHC, and pathways with Bayes factors greater than 1e10 (over
        % Bayes factor for enrichment of MHC) given enrichment of the MHC
        't1d' 't1d/pathway-t1d-J.mat' 1 2.70e54 [3:6 9 14:16]
        't1d' 't1d/pathway-t1d-G.mat' 0 1       [1 2 4]
        %
        % Top 3 pathways in RA without conditioning on enrichment of the
        % MHC, and pathways with Bayes factors greater than 100 (over
        % Bayes factor for enrichment of MHC) given enrichment of the MHC.
        'ra'  'ra/pathway-ra-J.mat'   1 2.43e21 1:4
        'ra'  'ra/pathway-ra-G.mat'   0 1      [1 2 5]
        %
        % Pathways with Bayes factors greater than 100 in CD, plus the
        % MHC and xMHC.
        'cd'  'cd/pathway-cd-D.mat'   0 1       [1 8 9 11 12 14 19:21]
        %
        % Pathways with Bayes factors greater than 100 in T2D.
        't2d' 't2d/pathway-t2d-E.mat' 0 1       1:5 
        %
        % Pathways with Bayes factors greater than 50 in CAD.
        'cad' 'cad/pathway-cad-E.mat' 0 1       1:7
        %
        % Pathways with Bayes factors greater than 10 in BD.
        'bd'  'bd/pathway-bd-D.mat'   0 1       1:4 };
        
  d             = D(:,1);
  resultsfile   = D(:,2);
  mhcisenriched = D(:,3);
  BFmhc         = cell2mat(D(:,4));
  hind          = D(:,5);

  % This is the pathway index corresponding to the major histocompatibility
  % complex (MHC).
  mhcind = 3322;

  figure(1);
  clf
  windowpos = get(gcf,'Position');
  set(gcf,'Name','Results (Part K)','Color','white',...
          'NumberTitle','off','PaperPositionMode','auto',...
          'Position',[windowpos(1:2) 1000 800]);

  % Get a couple pieces of information about the selected pathways.
  n = length(d);
  N = cellfun(@numel,hind);
  for i = 1:n
    d{i}             = repmat(d(i),N(i),1);
    mhcisenriched{i} = repmat(mhcisenriched{i},N(i),1);  
  end
  d             = vertcat(d{:});
  mhcisenriched = vertcat(mhcisenriched{:});

  % Initialize storage for the combined results.
  N            = sum(N);
  paths        = zeros(N,1);
  numgenes     = zeros(N,1);
  numsnps      = zeros(N,1);
  numgenesmhc  = zeros(N,1);
  numsnpsmhc   = zeros(N,1);
  E0           = zeros(N,1);
  a0           = zeros(N,1);
  b0           = zeros(N,1);
  E            = zeros(N,1);
  a            = zeros(N,1);
  b            = zeros(N,1);

  % Repeat for each set of results.
  % 
  %   i  set of results
  %   j  pathway
  %   k  selected enrichment hypothesis
  %   t  final entry in table or barplot
  % 
  t = 0;
  for i = 1:n

    % Load the Bayes factors and other statistics.
    load(strcat(resultsdir,'/',resultsfile{i}));

    % Select data about the requested enrichment hypotheses.
    k     = hind{i};
    BF    = BF(k);
    H     = H(:,k);
    logw1 = logw1(:,:,k);

    % Sort the enrichment hypotheses from largest Bayes factor to the
    % smallest. 
    [ans k] = sort(-BF);
    BF      = BF(k);
    H       = H(:,k);
    logw1   = logw1(:,:,k);

    % Get the number of selected enrichment hypotheses.
    m = length(BF);

    % Repeat for each enrichment hypothesis.
    for k = 1:m
      t = t + 1;

      % Get the index of the enriched pathway.
      j        = find(H(:,k));
      paths(t) = j;

    % Get the number of genes assigned to the enriched pathway, and the number
    % of genes assigned to the MHC if the MHC region is enriched for disease
    % associations. If the enrichment hypothesis includes enrichment of the
    % MHC, then remove all genes within the MHC from this count.
      if mhcisenriched(t)
        numgenesmhc(t) = sum(full(pathway.genes(:,mhcind) & inref(gene)));
        numgenes(t)    = sum(pathway.genes(:,j) & ...
                             ~pathway.genes(:,mhcind) & ...
                             inref(gene));
      else
        numgenesmhc(t) = 0;
        numgenes(t)    = sum(pathway.genes(:,j) & inref(gene));
      end

      % Get the number of SNPs assigned to the enriched pathway, and the
      % number of SNPs assigned to the MHC if the MHC region is enriched
      % for disease associations. If the enriched hypothesis includes
      % enrichment of the MHC, then remove all SNPs within the MHC from
      % this count.
      if mhcisenriched(t)
        numsnpsmhc(t) = sum(A(:,mhcind));
        numsnps(t)    = sum(A(:,j) & ~A(:,mhcind));
      else
        numsnpsmhc(t) = 0;
        numsnps(t)    = sum(A(:,j));
      end

      % Get the posterior mean and 95% credible interval interval for the
      % genome-wide log-odds (theta0).
      w             = normalizelogweights(logw1(:,:,k));
      E0(t)         = dot(theta0,sumcols(w));
      [a0(t) b0(t)] = cred(theta0,sumcols(w),E0(t));
      
      % Get the posterior mean and 95% credible interval interval for the
      % enrichment parameter (theta)
      E(t)        = dot(theta,sumrows(w));
      [a(t) b(t)] = cred(theta,sumrows(w),E(t));

      % Print more details about the enrichment hypothesis.
      if isempty(pathway.synonyms{j})
        star = '';
      else
        star = '*';
      end
      if mhcisenriched(t)
        fprintf('%s: %s + MHC:\n',upper(d{t}),pathway.nickname{j});
      else
        fprintf('%s: %s\n',upper(d{t}),pathway.nickname{j});
      end
      fprintf('%s%s (%s/%s,i=%d)\n',pathway.label{j},star,...
              pathway.source{j},pathway.database{j},j);
      bf = BF(k) * BFmhc(i);
      if bf > 1e4
        fprintf('BF     = %0.2e\n',bf);
      else
        fprintf('BF     = %d\n',round(bf));
      end
      fprintf('theta0 = %0.2f (%0.1f,%0.1f)\n',a0(t),E0(t),b0(t));
      fprintf('theta  = %0.2f (%0.1f,%0.1f)\n',a(t),E(t),b(t));
      fprintf('\n');
    end
  end

  % Plot the number of genes.
  subplot(1,6,2:3);
  h = barh(1:N,numgenesmhc + numgenes,0.8);
  set(h,'EdgeColor',rgb('white'),'LineWidth',3,'FaceColor',rgb('darkorange'));
  hold on
  h = barh(1:N,numgenesmhc,0.8);
  set(h,'EdgeColor',rgb('white'),'LineWidth',3,'FaceColor',rgb('darkorange'));
  hold off
  set(gca,'FontSize',9,'FontName','fixed');
  set(gca,'YDir','reverse','YLim',[0 N+1],'YTick',1:N);
  set(gca,'XLim',[0 800],'XTick',0:200:800);
  set(gca,'TickLength',[0 0]);
  xlabel('no. genes','FontSize',10);
  str = num2str(numgenes,'%-d');
  I   = find(~mhcisenriched);
  text(numgenes(I) + 10,I,str(I,:),'FontSize',9,'FontName','fixed');
  str = strcat(num2str(numgenesmhc),'+',str);
  I   = find(mhcisenriched);
  text(numgenesmhc(I) + numgenes(I) + 10,I,str(I,:),...
       'FontSize',9,'FontName','fixed');

  % Label this table entry by the enriched pathway(s).
  str    = strcat(upper(d),'-',pathway.nickname(paths));
  I      = find(mhcisenriched);
  str(I) = strcat(str(I),'+','MHC');
  set(gca,'YTickLabel',str);

  % Plot the number of SNPs.
  subplot(1,6,4);
  h = barh(1:N,numsnpsmhc + numsnps,0.7);
  set(h,'LineStyle','none','FaceColor',rgb('darkorange'));
  hold on
  h = barh(1:N,numsnpsmhc,0.7);
  set(h,'LineStyle','none','FaceColor',rgb('cornflowerblue'));
  hold off
  set(gca,'FontSize',9,'FontName','fixed');
  set(gca,'YDir','reverse','YLim',[0 N+1],'YTick',[]);
  set(gca,'XLim',[0 2.5e4],'XTick',[0:1e4:2e4]);
  xlabel('no. SNPs','FontSize',10);
  str = num2str(numsnps,'%-d');
  I   = find(~mhcisenriched);
  text(numsnps(I) + 1e3,I,str(I,:),'FontSize',9,'FontName','fixed');
  str = strcat(num2str(numsnpsmhc),'+',str);
  I   = find(mhcisenriched);
  text(numsnpsmhc(I) + numsnps(I) + 1e3,I,str(I,:),...
       'FontSize',9,'FontName','fixed');

  % Show the posterior mean and 95% credible interval of the genome-wide
  % log-odds (theta0) for each enrichment hypothesis.
  subplot(1,6,5);
  h = herrorbar(E0,1:N,E0-a0,b0-E0,'ko',0.3);
  set(h,'Markersize',4,'MarkerFaceColor',rgb('black'));
  set(gca,'FontSize',9,'FontName','fixed');
  set(gca,'YDir','reverse','YLim',[0 N+1],'YTick',[]);
  set(gca,'XLim',[-6.25 -2],'XTick',-6:-2,'XGrid','on');
  set(gca,'TickLength',[0 0]);
  xlabel('genome-wide log-odds','FontSize',10);

  % Show the posterior mean and 95% credible interval of the enrichment
  % parameter (theta) for each enrichment hypothesis.
  subplot(1,6,6);
  h = herrorbar(E,1:N,E-a,b-E,'ko',0.3);
  set(h,'Markersize',4,'MarkerFaceColor',rgb('black'));
  set(gca,'FontSize',9,'FontName','fixed');
  set(gca,'YDir','reverse','YLim',[0 N+1],'YTick',[]);
  set(gca,'XLim',[0 5],'XTick',0:5,'XGrid','on');
  set(gca,'TickLength',[0 0]);
  xlabel('log_{10}enrichment','FontSize',10);

 case 'L'

  % SCATTERPLOT OF CD ASSOCIATIONS FOR CYTOKINE SIGNALING
  % -----------------------------------------------------
  % These are labels corresponding to the selected segments.
  seglabel = cell(2e4,1);
  seglabel{403}   = 'IL23R';          
  seglabel{2873}  = 'ATG16L1(SHIP1)';
  seglabel{5599}  = 'PTGER4';         
  seglabel{6125}  = 'IBD5(IRF1)';           
  seglabel{6320}  = 'IL12B';
  seglabel{6718}  = 'MHC(class I)';  
  seglabel{6745}  = 'MHC(class II)'; 
  seglabel{10868} = 'ZNF365';         
  seglabel{10945} = 'CAMK2G';         
  seglabel{11134} = 'NKX2-3';         
  seglabel{15256} = 'NOD2';           
  seglabel{15545} = 'IRF8';           
  seglabel{15764} = 'STAT3';          
  seglabel{16078} = 'PTPN2';          

  figure(1);
  clf
  windowpos = get(gcf,'Position');
  set(gcf,'Name','Results (Part L)','Color','white',...
          'NumberTitle','off','PaperPositionMode','auto',...
          'Position',[windowpos(1:2) 640 570]);

  % Load the genotype and phenotype data.
  fprintf('Loading genotype and phenotype data for CD.\n');
  load(strcat(datadir,'/cd.mat'));

  % Create overlapping segments with 50 SNPs.
  fprintf('Creating segments.\n');
  [Aseg segchr segpos] = makesegments2(chr,pos,25,5e4);
  [Aseg segchr segpos] = combinesegments(Aseg,segchr,segpos);

  % NULL HYPOTHESIS - NO ENRICHED PATHWAYS
  % Load the results conditioned on the null hypothesis that no pathways are
  % enriched and, for each segment, compute the posterior statistics P1 and
  % P2 averaged over settings of the hyperparameters.
  fprintf('NULL HYPOTHESIS - NO ENRICHED PATHWAYS.\n');
  fprintf('Loading null hypothesis results.\n');
  null            = load(strcat(resultsdir,'/cd/pathway-cd-B.mat'));
  w0              = normalizelogweights(null.logw);
  [P1null P2null] = compilesegstats(Aseg,null.alpha,w0);

  % ENRICHMENT OF CYTOKINE SIGNALING GENES
  % Load the results conditioned on the hypothesis that cytokine signaling
  % genes are enriched.
  fprintf('ENRICHMENT OF CYTOKINE SIGNALING GENES.\n');
  fprintf('Loading pathway enrichment results.\n');
  load(strcat(resultsdir,'/cd/pathway-cd-K.mat'));
  i         = 3;
  n1        = numel(theta);
  logw1     = logw1(:,:,i);
  alphapath = alpha{i};
  mupath    = mu{i};

  % Re-load the results conditioned on the null hypothesis that no pathways
  % are enriched, this time selecting the hyperparameter settings that
  % match the hyperparameter settings under the alternative hypothesis.
  null = getnullstats(strcat(resultsdir,'/cd/pathway-cd-B.mat'),theta0);

  % Get the index of the enriched pathway (j), and the set of SNPs assigned to
  % the enriched pathway (snps). Then get the posterior inclusion
  % probabilities (and posterior mean additive effects) for each setting of
  % the hyperparameters.
  j     = find(H(:,i));
  snps  = find(A(:,j));
  alpha = repmat(null.alpha,1,n1);
  mu    = repmat(null.mu,1,n1);

  alpha(snps,:) = alphapath(snps,:);
  mu(snps,:)    = mupath(snps,:);

  % Compute the posterior inclusion probabilities (PIPs) and posterior mean
  % regression coefficients averaged over the settings of the
  % hyperparameters and, for each segment, compute the posterior statistics
  % P1 and P2 averaged over settings of the hyperparameters.
  w       = normalizelogweights(logw1);
  PIP     = alpha * w(:);
  mu      = mu * w(:);
  [P1 P2] = compilesegstats(Aseg,alpha,w(:));  

  % Plot the posterior probabilities P1 given that no pathways are enriched
  % (horizontal axis) versus the posterior probabilities P1 given that
  % cytokine signaling genes are enriched (vertical axis). I highlight
  % segments that overlap SNPs assigned to the enriched pathway.
  plot([0 2],[0 2],':','Color',rgb('dodgerblue'),'LineWidth',1);
  hold on
  inpath = full(sumrows(Aseg(snps,:)))' > 0;

  % First plot "non-redundant" segments that have posterior probabilities
  % less than 0.4 under the alternative hypothesis.
  I = (P1 > 0.01 | P1null > 0.01) & P1 < 0.4 & ~redundantsegment(Aseg,P1);
  t = I & ~inpath;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('black'),...
       'MarkerEdgeColor','none','MarkerSize',6);
  t = I & inpath;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('firebrick'),...
       'MarkerEdgeColor','none','MarkerSize',6);

  % Next plot "non-redundant" segment with posterior probabilities greater
  % than 0.4 under the alternative hypothesis.
  I = P1 >= 0.4 & ~redundantsegment(Aseg,P1);
  t = I & ~inpath;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('black'),...
       'MarkerEdgeColor','none','MarkerSize',6);
  t = I & inpath;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('firebrick'),...
       'MarkerEdgeColor','none','MarkerSize',6);

  % Label some of the segments.
  t = I & ~inpath;
  text(P1null(t) + 0.005,P1(t) + 0.005,seglabel(t,:),...
       'VerticalAlignment','bottom','FontSize',9,...
       'FontName','fixed','Color',rgb('black'));
  t = I & inpath;
  text(P1null(t) + 0.005,P1(t) + 0.005,seglabel(t,:),...
       'VerticalAlignment','bottom','FontSize',9,...
       'FontName','fixed','Color',rgb('firebrick'));

  hold off
  set(gca,'FontSize',10,'FontName','fixed');
  set(gca,'XLim',[0 1.05],'XTick',0:0.5:1);
  set(gca,'YLim',[0 1.05],'YTick',0:0.5:1);
  set(gca,'TickDir','out');
  xlabel('P1 when no pathways are enriched');
  ylabel('P1 when cytokine signaling genes are enriched');
  title('Crohn''s disease');

  % SCATTERPLOT OF CD ASSOCIATIONS FOR 2 ENRICHED PATHWAYS
  % ------------------------------------------------------
  figure(2);
  clf
  windowpos = get(gcf,'Position');
  set(gcf,'Name','Results (Part L-2)','Color','white',...
          'NumberTitle','off','PaperPositionMode','auto',...
          'Position',[windowpos(1:2) 1050 425]);

  % Save the results for enrichment of cytokine signaling.
  P1null = P1;
  P2null = P2;

  % Load the results conditioned on enrichment of 2 pathways.
  fprintf('ENRICHMENT OF 2 PATHWAYS.\n');
  fprintf('Loading pathway enrichment results.\n');
  load(strcat(resultsdir,'/cd/pathway-cd-G.mat'));
  n1        = numel(theta);
  alphapath = alpha;

  % Re-load the results conditioned on the null hypothesis that no pathways
  % are enriched, this time selecting the hyperparameter settings that match
  % the hyperparameter settings under the alternative hypothesis.
  null = getnullstats(strcat(resultsdir,'/cd/pathway-cd-B.mat'),theta0);

  % Compute the posterior statistics P1 for each enrichment hypothesis.
  n    = length(BF);
  nseg = length(segchr);
  P1   = zeros(nseg,n);
  P2   = zeros(nseg,n);
  fprintf('Computing segment statistics given enrichment hypotheses:\n');
  for i = 1:n
    
    % Get the indices of the enriched pathways (paths), and the set of SNPs
    % assigned to the enriched pathways (snps).
    paths = find(H(:,i))';
    snps  = find(sumcols(A(:,paths)));
    fprintf('%d. Enrichment of "%s" + "%s"\n',i,...
            pathway.label{paths(1)},...
            pathway.label{paths(2)});

    % Get the posterior inclusion probabilities for each setting of the
    % hyperparameters.
    alpha         = repmat(null.alpha,[1 n1]);
    alpha(snps,:) = alphapath{i}(snps,:);

    % For each segment, compute the posterior statistics P1 and P2
    % averaged over settings of the hyperparameters.
    w                 = normalizelogweights(logw1(:,:,i));
    [P1(:,i) P2(:,i)] = compilesegstats(Aseg,alpha,w(:));
  end

  % Get the final posterior statistics P1 and P2 by averaging over the
  % enrichment hypotheses.
  P1 = P1 * BF / sum(BF);
  P2 = P2 * BF / sum(BF);

  % Plot the posterior probabilities given that cytokine signaling genes are
  % enriched (horizontal axis) versus the posterior probabilities given that
  % 2 pathways are enriched (vertical axis). I highlight segments that
  % overlap SNPs assigned to the additional enriched pathways (and are not
  % assigned to the "cytokine signaling in immune system" patwhay).
  j      = 2069;
  paths  = find(sumcols(H));
  snps   = find(sumcols(A(:,paths)) & ~A(:,j));
  inpath = full(sumrows(Aseg(snps,:)))' > 0;
  subplot(1,2,1);
  plot([0 2],[0 2],':','Color',rgb('dodgerblue'),'LineWidth',1);
  hold on

  % Plot "non-redundant" segments.
  I = (P1 > 0.01 | P1null > 0.01) & ~redundantsegment(Aseg,P1);
  t = I & ~inpath;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('black'),...
       'MarkerEdgeColor','none','MarkerSize',6);
  t = I & inpath;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('firebrick'),...
       'MarkerEdgeColor','none','MarkerSize',6);

  % Label some of the segments.
  t = I & inpath;
  text(P1null(t) + 0.005,P1(t) + 0.005,seglabel(t,:),...
       'VerticalAlignment','bottom','FontSize',9,...
       'FontName','fixed','Color',rgb('firebrick'));

  hold off
  set(gca,'FontSize',10,'FontName','fixed');
  set(gca,'XLim',[0 1.05],'XTick',0:0.5:1);
  set(gca,'YLim',[0 1.05],'YTick',0:0.5:1);
  set(gca,'TickDir','out');
  xlabel('P1 when cytokine signaling genes are enriched');
  ylabel('P1 when 2 pathways are enriched');
  title('Crohn''s disease');

  % HYPOTHESES IN WHICH 3 PATHWAYS ARE ENRICHED
  % -------------------------------------------
  % Save the results for enrichment of 2 pathways.
  P1null = P1;
  P2null = P2;

  % Load the results conditioned on enrichment of 3 pathways.
  fprintf('ENRICHMENT OF 3 PATHWAYS.\n');
  fprintf('Loading pathway enrichment results.\n');
  load(strcat(resultsdir,'/cd/pathway-cd-J.mat'));
  n1        = numel(theta);
  alphapath = alpha;

  % Re-load the results conditioned on the null hypothesis that no pathways
  % are enriched, this time selecting the hyperparameter settings that match
  % the hyperparameter settings under the alternative hypothesis.
  null = getnullstats(strcat(resultsdir,'/cd/pathway-cd-B.mat'),theta0);

  % Compute the posterior statistics P1 for each enrichment hypothesis.
  n    = length(BF);
  nseg = length(segchr);
  P1   = zeros(nseg,n);
  P2   = zeros(nseg,n);
  fprintf('Computing segment statistics given enrichment hypotheses:\n');
  for i = 1:n

    % Get the index of the enriched pathways (paths), and the set of SNPs
    % assigned to the enriched pathways (snps).
    paths = find(H(:,i))';
    snps  = find(sumcols(A(:,paths)));
    fprintf('%d. Enrichment of "%s" + "%s" + "%s"\n',i,...
            pathway.label{paths(1)},...
            pathway.label{paths(2)},...
            pathway.label{paths(3)});
    
    % Get the posterior inclusion probabilities for each setting of the
    % hyperparameters.
    alpha         = repmat(null.alpha,[1 n1]);
    alpha(snps,:) = alphapath{i}(snps,:);

    % For each segment, compute the posterior statistics P1 and P2
    % averaged over settings of the hyperparameters.
    w                 = normalizelogweights(logw1(:,:,i));
    [P1(:,i) P2(:,i)] = compilesegstats(Aseg,alpha,w(:));
  end

  % Get the final posterior statistics P1 and P2 by averaging over the
  % enrichment hypotheses.
  P1 = P1 * BF / sum(BF);
  P2 = P2 * BF / sum(BF);

  % Plot the posterior probabilities given that 2 pathways are enriched
  % (horizontal axis) versus the posterior probabilities given that 3
  % pathways are enriched (vertical axis). I highlight segments that overlap
  % SNPs assigned to the additional enriched pathways (and are not assigned
  % to either the cytokine signaling or the IL-23 signaling pathways).
  paths  = find(sumcols(H));
  snps   = find(sumcols(A(:,paths)) & ~A(:,2069) & A(:,2591));
  inpath = full(sumrows(Aseg(snps,:)))' > 0;
  subplot(1,2,2);
  plot([0 2],[0 2],':','Color',rgb('dodgerblue'),'LineWidth',1);
  hold on

  % Plot "non-redundant" segments.
  I = (P1 > 0.01 | P1null > 0.01) & ~redundantsegment(Aseg,P1);
  t = I & ~inpath;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('black'),...
       'MarkerEdgeColor','none','MarkerSize',6);
  t = I & inpath;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('firebrick'),...
       'MarkerEdgeColor','none','MarkerSize',6);

  % Label some of the segments.
  t = I & inpath;
  text(P1null(t) + 0.005,P1(t) + 0.005,seglabel(t,:),...
       'VerticalAlignment','bottom','FontSize',9,...
       'FontName','fixed','Color',rgb('firebrick'));

  hold off
  set(gca,'FontSize',10,'FontName','fixed');
  set(gca,'XLim',[0 1.05],'XTick',0:0.5:1);
  set(gca,'YLim',[0 1.05],'YTick',0:0.5:1);
  set(gca,'TickDir','out');
  xlabel('P1 when 2 pathways are enriched');
  ylabel('P1 when 3 pathways are enriched');
  title('Crohn''s disease');

 case 'M'

  % SELECTED REGIONS OF GENOME FOR CD GIVEN ENRICHMENT OF CYTOKINE SIGNALING
  % ------------------------------------------------------------------------
  % The columns of this table specify: (1-3) the chromosome and the start and
  % end positions on the chromosome (in Mb) that make up the "critical
  % region" for the disease association; (4) candidate gene(s) corresponding
  % to locus.
  D = {  1  67.30   67.48 'IL23R'
         2 233.92  234.27 'ATG16L1'
         5  40.32   40.66 'PTGER4'
         5 129.54  132.04 'IBD5'
         6  25.52   33.76 'MHC'
        10  64.00   64.43 'ZNF365'
        10  101.26 101.32 'NKX2-3'
        16  49.00   49.40 'NOD2'
        17  37.50   38.30 'STAT3'
        18  12.76   12.91 'PTPN2' };

  region.chr   = cell2mat(D(:,1));
  region.start = cell2mat(D(:,2)) * 1e6;
  region.stop  = cell2mat(D(:,3)) * 1e6;
  region.label = D(:,4);

  % Load the genotype and phenotype data.
  fprintf('Loading genotype and phenotype data for CD.\n');
  load(strcat(datadir,'/cd.mat'));

  % Get the cases and controls.
  cases = find(y == 1);
  ctrls = find(y == 0);

  % Get the set of SNPs for each selected region.
  Ar = snps2genes(chr,pos,region,0,1e4);

  % NULL HYPOTHESIS - NO ENRICHED PATHWAYS
  % Load the results conditioned on the null hypothesis that no pathways are
  % enriched and compute: (1) for each of the selected regions, the
  % posterior statistics P1 and P2 averaged over settings of the
  % hyperparameters; (2) for each SNP, the posterior inclusion probability
  % averaged over settings of the hyperparameters.
  fprintf('Loading null hypothesis results.\n');
  null            = load(strcat(resultsdir,'/cd/pathway-cd-B.mat'));
  w0              = normalizelogweights(null.logw);
  PIP0            = null.alpha * w0;
  [P1null P2null] = compilesegstats(Ar,null.alpha,w0);

  % Load the results conditioned on the hypothesis that cytokine signaling
  % genes are enriched.
  fprintf('Loading pathway enrichment results.\n');
  load(strcat(resultsdir,'/cd/pathway-cd-K.mat'));
  i         = 3;
  n1        = numel(theta);
  logw1     = logw1(:,:,i);
  alphapath = alpha{i};
  mupath    = mu{i};
  spath     = s{i};

  % ENRICHMENT OF CYTOKINE SIGNALING GENES
  % Re-load the results conditioned on the null hypothesis that no pathways
  % are enriched, this time selecting the hyperparameter settings that match
  % the hyperparameter settings under the alternative hypothesis.
  null = getnullstats(strcat(resultsdir,'/cd/pathway-cd-B.mat'),theta0);

  % Get the index of the enriched pathway (j), and the set of SNPs assigned to
  % the enriched pathway (snps). Then get the posterior inclusion
  % probabilities and other posterior expectations for each setting of the
  % hyperparameters.
  j     = find(H(:,i));
  snps  = find(A(:,j));
  alpha = repmat(null.alpha,1,n1);
  mu    = repmat(null.mu,1,n1);
  s     = repmat(null.s,1,n1);

  alpha(snps,:) = alphapath(snps,:);
  mu(snps,:)    = mupath(snps,:);
  s(snps,:)     = spath(snps,:);

  % Compute the posterior inclusion probabilities (PIPs) averaged over the
  % settings of the hyperparameters and, for each selected region, compute
  % the posterior statistics P1 and P2 averaged over settings of the
  % hyperparameters.
  w       = normalizelogweights(logw1);
  PIP     = alpha * w(:);
  [P1 P2] = compilesegstats(Ar,alpha,w(:));  

  % Show the table legend.
  fprintf('RESULTS (PART M)\n');
  fprintf('                        P1   P1   P2   P2         ');
  fprintf('                                  MAF   MAF\n');
  fprintf('dis chr  region (Mb)  null  alt null  alt gene(s) ');
  fprintf('SNP         PIP   LOR(95%% CI)    ctrls cases\n');

  % Repeat for each selected region.
  n = length(region.label);
  for i = 1:n

    % Get the set of SNPs in the region.
    snps = find(Ar(:,i));
          
    % Get the SNP in the region with the largest posterior inclusion
    % probability.
    [ans k] = max(PIP(snps));
    k       = snps(k);
      
    % Get the mean and 95% credible interval for the additive effect
    % (i.e. log-odds ratio) corresponding to the SNP with the largest
    % PIP in the segment.
    E     = dot(w(:),mu(k,:));
    beta  = (-2:0.01:2)';
    r     = normpdf(repmat(beta,1,numel(w)),...
                    repmat(mu(k,:),length(beta),1),...
                    repmat(s(k,:),length(beta),1));
    [a b] = cred(beta,r*w(:),E);

    % Print statistics about the selected region.
    fprintf(' CD %3d %6.2f-%-6.2f %0.2f %0.2f %0.2f %0.2f %7s ',...
            region.chr(i),region.start(i)/1e6,region.stop(i)/1e6,...
            P1null(i),P1(i),P2null(i),P2(i),region.label{i});
    fprintf('%-10s %0.2f ',sprintf('rs%d',labels(k)),PIP(k));
    fprintf('%+0.2f(%0.2f-%0.2f) %0.3f %0.3f\n',...
            E,min(abs([a b])),max(abs([a b])),...
            maf(X(ctrls,k)),maf(X(cases,k)));
  end

 case 'N'

  % SCATTERPLOT OF RA ASSOCIATIONS FOR "MEASLES"
  % --------------------------------------------
  % These are labels corresponding to selected segments.
  seglabel = cell(2e4,1);
  seglabel{7}     = 'TP73';
  seglabel{702}   = 'PTPN22';
  seglabel{2605}  = 'STAT4';
  seglabel{3931}  = 'IL12A';
  seglabel{6305}  = 'IL12B';
  seglabel{9584}  = 'JAK2';
  seglabel{9735}  = 'IFNB1';
  seglabel{7367}  = 'TNFAIP3';
  seglabel{8362}  = 'IRF5';
  seglabel{10420} = 'IL2RA';
  seglabel{10424} = 'PRKCQ';
  seglabel{11685} = 'TRAF6';
  seglabel{12711} = 'KIF5A(CDK4)';
  seglabel{16510} = 'TYK2,PIN1';
  seglabel{17305} = 'IFNAR1';
  seglabel{17550} = 'IL2RB';

  figure(1);
  clf
  windowpos = get(gcf,'Position');
  set(gcf,'Name','Results (Part N)','Color','white',...
          'NumberTitle','off','PaperPositionMode','auto',...
          'Position',[windowpos(1:2) 1050 425]);

  % Load the genotype and phenotype data.
  fprintf('Loading genotype and phenotype data for RA.\n');
  load(strcat(datadir,'/ra.mat'));

  % Create overlapping segments with 50 SNPs.
  fprintf('Creating segments.\n');
  [Aseg segchr segpos] = makesegments2(chr,pos,25,5e4);
  [Aseg segchr segpos] = combinesegments(Aseg,segchr,segpos);

  % NULL HYPOTHESIS - NO ENRICHED PATHWAYS
  % Load the results conditioned on the null hypothesis that no pathways are
  % enriched and, for each segment, compute the posterior statistics P1 and
  % P2 averaged over settings of the hyperparameters.
  fprintf('Loading null hypothesis results.\n');
  null            = load(strcat(resultsdir,'/ra/pathway-ra-D.mat'));
  w0              = normalizelogweights(null.logw);
  [P1null P2null] = compilesegstats(Aseg,null.alpha,w0);

  % ENRICHMENT OF "MEASLES" + MHC
  % Load the results conditioned on the hypothesis that the Measles pathway
  % and SNPs within the MHC region are enriched.
  fprintf('Loading pathway enrichment results.\n');
  load(strcat(resultsdir,'/ra/pathway-ra-J.mat'));
  i         = 3;
  n1        = numel(theta);
  logw1     = logw1(:,:,i);
  alphapath = alpha{i};
  mupath    = mu{i};

  % Re-load the results conditioned on the null hypothesis that the MHC is
  % enriched, this time selecting the hyperparameter settings that match the
  % hyperparameter settings under the alternative hypothesis.
  null = getnullstats(strcat(resultsdir,'/ra/pathway-ra-D.mat'),theta0);
  mhcnull = load(strcat(resultsdir,'/ra/pathway-ra-G.mat'));
  [logw0 alpha0 mu0 s0] = getpathwaynull(null,mhcnull,theta0,thetaMHC,1);

  % Get the SNPs that lie within the MHC region.
  p = length(labels);
  isinmhc = zeros(p,1);
  isinmhc(mhcsnps) = 1;

  % Get the index of the enriched pathway (j), and the set of SNPs assigned to
  % the enriched pathway that are not in the MHC (snps). Then get the
  % posterior inclusion probabilities (and posterior mean additive effects)
  % for each setting of the hyperparameters.
  j             = find(H(:,i));
  snps          = find(A(:,j) & ~isinmhc);
  alpha         = repmat(alpha0,1,n1);
  mu            = repmat(mu0,1,n1);
  alpha(snps,:) = alphapath(snps,:);
  mu(snps,:)    = mupath(snps,:);

  % Compute the posterior inclusion probabilities (PIPs) and posterior mean
  % regression coefficients averaged over the settings of the
  % hyperparameters and, for each segment, compute the posterior statistics
  % P1 and P2 averaged over settings of the hyperparameters.
  w       = normalizelogweights(logw1);
  PIP     = alpha * w(:);
  mu      = mu * w(:);
  [P1 P2] = compilesegstats(Aseg,alpha,w(:));

  % Plot the posterior probabilities given that no pathways are enriched
  % (horizontal axis) versus the posterior probabilities given that the
  % "Measles" pathway, including the MHC, are enriched (vertical axis). I
  % highlight segments that overlap SNPs assigned to the enriched pathway.
  subplot(1,2,1);
  plot([0 2],[0 2],':','Color',rgb('dodgerblue'),'LineWidth',1);
  hold on
  inpathseg = full(sumrows(Aseg(snps,:)))' > 0;
  inmhcseg  = full(sumrows(Aseg(mhcsnps,:)))' > 0;

  % First, plot "non-redundant" segments overlapping the MHC.
  I = (P1 > 0.01 | P1null > 0.01) & inmhcseg & ~redundantsegment(Aseg,P1);
  plot(P1null(I),P1(I),'o','MarkerFaceColor',rgb('royalblue'),...
       'MarkerEdgeColor','none','MarkerSize',6);

  % Next, plot "non-redundant" segments that have posterior probabilities
  % less than 0.4 under the alternative hypothesis.
  I = (P1 > 0.01 | P1null > 0.01) & P1 < 0.4 & ~inmhcseg & ...
      ~redundantsegment(Aseg,P1);
  t = I & ~inpathseg;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('black'),...
       'MarkerEdgeColor','none','MarkerSize',6);
  t = I & inpathseg;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('firebrick'),...
       'MarkerEdgeColor','none','MarkerSize',6);

  % Finally, plot "non-redundant" segments with posterior probabilities
  % greater than 0.4 under the alternative hypothesis. Here I need to fix
  % one segment at 6q23 that causes a bit of a problem because it
  % underestimates the posterior probability under the null hypothesis.
  I       = P1 >= 0.4 & ~inmhcseg & ~redundantsegment(Aseg,P1);
  I(7367) = 1;
  I(7368) = 0;

  t = I & ~inpathseg;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('black'),...
       'MarkerEdgeColor','none','MarkerSize',6);
  t = I & inpathseg;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('firebrick'),...
       'MarkerEdgeColor','none','MarkerSize',6);

  % Label some of the segments.
  I = (P1 > 0.01 | P1null > 0.01) & ~inmhcseg & ...
      ~redundantsegment(Aseg,P1);
  t = I & ~inpathseg;
  text(P1null(t) + 0.005,P1(t) + 0.005,seglabel(t,:),...
       'VerticalAlignment','bottom','FontSize',9,'FontName','fixed',...
       'Color',rgb('black'));
  t = I & inpathseg;
  text(P1null(t) + 0.005,P1(t) + 0.005,seglabel(t,:),...
       'VerticalAlignment','bottom','FontSize',9,'FontName','fixed',...
       'Color',rgb('firebrick'));

  hold off
  set(gca,'FontSize',10,'FontName','fixed');
  set(gca,'XLim',[0 1.05],'XTick',0:0.5:1);
  set(gca,'YLim',[0 1.05],'YTick',0:0.5:1);
  set(gca,'TickDir','out');
  xlabel('P1 when no pathways are enriched');
  ylabel('P1 when MHC and "Measles" pathway are enriched');
  title('Rheumatoid arthritis');

  % SCATTERPLOT OF RA ASSOCIATIONS FOR 2 ENRICHED PATHWAYS
  % ------------------------------------------------------
  % Save the results for enrichment of "Measles" + MHC.
  P1null = P1;
  P2null = P2;

  % Load the results conditioned on enrichment of 2 pathways.
  fprintf('ENRICHMENT OF 2 PATHWAYS.\n');
  fprintf('Loading pathway enrichment results.\n');
  load(strcat(resultsdir,'/ra/pathway-ra-L.mat'));
  n1        = numel(theta);
  alphapath = alpha;

  % Re-load the results conditioned on the null hypothesis that the MHC is
  % enriched, this time selecting the hyperparameter settings that match the
  % hyperparameter settings under the alternative hypothesis.
  null = getnullstats(strcat(resultsdir,'/ra/pathway-ra-D.mat'),theta0);
  mhcnull = load(strcat(resultsdir,'/ra/pathway-ra-G.mat'));
  [logw0 alpha0 mu0 s0] = getpathwaynull(null,mhcnull,theta0,thetaMHC,1);

  % Compute the posterior statistics P1 for each enrichment hypothesis.
  n    = length(BF);
  nseg = length(segchr);
  P1   = zeros(nseg,n);
  P2   = zeros(nseg,n);
  fprintf('Computing segment statistics given enrichment hypotheses:\n');
  for i = 1:n

    % Get the indices of the enriched pathways (paths), and the set of SNPs
    % assigned to the enriched pathways (snps).
    paths = find(H(:,i))';
    snps  = find(sumcols(A(:,paths)));
    fprintf('%d. Enrichment of MHC + "%s" + "%s"\n',i,...
            pathway.label{paths(1)},...
            pathway.label{paths(2)});

    % Get the posterior inclusion probabilities for each setting of the
    % hyperparameters.
    alpha         = repmat(alpha0,[1 n1]);
    alpha(snps,:) = alphapath{i}(snps,:);

    % For each segment, compute the posterior statistics P1 and P2
    % averaged over settings of the hyperparameters.
    w                 = normalizelogweights(logw1(:,:,i));
    [P1(:,i) P2(:,i)] = compilesegstats(Aseg,alpha,w(:));
  end

  % Get the final posterior statistics P1 and P2 by averaging over the
  % enrichment hypotheses.
  P1 = P1 * BF / sum(BF);
  P2 = P2 * BF / sum(BF);

  % Plot the posterior probabilities given that "Measles" + MHC are enriched
  % (horizontal axis) versus the posterior probabilities given that 2
  % pathways are enriched (vertical axis). I highlight segments that overlap
  % SNPs assigned to the additional enriched pathways (and are not assigned
  % to the "Measles" pathway, nor to the MHC).
  j         = 2088;
  paths     = find(sumcols(H));
  snps      = find(sumcols(A(:,paths)) & ~A(:,j) & ~isinmhc);
  inpathseg = full(sumrows(Aseg(snps,:)))' > 0;
  inmhcseg  = full(sumrows(Aseg(mhcsnps,:)))' > 0;
  subplot(1,2,2);
  plot([0 2],[0 2],':','Color',rgb('dodgerblue'),'LineWidth',1);
  hold on

  % First, plot "non-redundant" segments overlapping the MHC.
  I = (P1 > 0.01 | P1null > 0.01) & inmhcseg & ~redundantsegment(Aseg,P1);
  plot(P1null(I),P1(I),'o','MarkerFaceColor',rgb('royalblue'),...
       'MarkerEdgeColor','none','MarkerSize',6);

  % Plot "non-redundant" segments that do not overlap the MHC.
  I = (P1 > 0.01 | P1null > 0.01) & ~inmhcseg & ~redundantsegment(Aseg,P1);
  t = I & ~inpathseg;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('black'),...
       'MarkerEdgeColor','none','MarkerSize',6);
  t = I & inpathseg;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('firebrick'),...
       'MarkerEdgeColor','none','MarkerSize',6);

  % Label some of the segments.
  t = I & inpathseg;
  text(P1null(t) + 0.005,P1(t) + 0.005,seglabel(t,:),...
       'VerticalAlignment','bottom','FontSize',9,...
       'FontName','fixed','Color',rgb('firebrick'));

  hold off
  set(gca,'FontSize',10,'FontName','fixed');
  set(gca,'XLim',[0 1.05],'XTick',0:0.5:1);
  set(gca,'YLim',[0 1.05],'YTick',0:0.5:1);
  set(gca,'TickDir','out');
  xlabel('P1 when "Measles" genes + MHC are enriched');
  ylabel('P1 when 2 pathways + MHC are enriched');
  title('rheumatoid arthritis');

 case 'O'

  % SELECTED REGIONS OF GENOME FOR RA GIVEN ENRICHMENT OF "MEASLES"
  % ---------------------------------------------------------------
  % The columns of this table specify: (1-3) the chromosome and the start and
  % end positions on the chromosome (in Mb) that make up the "critical
  % region" for the disease association; (4) candidate gene(s) corresponding
  % to locus.
  D = {  1   3.50   3.70 'TP73' 
         1 113.53 114.36 'PTPN22' 
         6  25.52  33.76 'MHC'
         6 138.00 138.47 'TNFAIP3'
        10   6.07   6.26 'IL2RA'
        10   6.36   6.49 'PRKCQ'
        12  55.77  56.82 'KIF5A'
        22  35.57  35.90 'IL2RB' };

  region.chr   = cell2mat(D(:,1));
  region.start = cell2mat(D(:,2)) * 1e6;
  region.stop  = cell2mat(D(:,3)) * 1e6;
  region.label = D(:,4);

  % Load the genotype and phenotype data.
  fprintf('Loading genotype and phenotype data for RA.\n');
  load(strcat(datadir,'/ra.mat'));

  % Get the cases and controls.
  cases = find(y == 1);
  ctrls = find(y == 0);

  % Get the set of SNPs for each selected region.
  Ar = snps2genes(chr,pos,region,0,1e4);

  % NULL HYPOTHESIS - NO ENRICHED PATHWAYS
  % Load the results conditioned on the null hypothesis that no pathways are
  % enriched and compute: (1) for each of the selected regions, the
  % posterior statistics P1 and P2 averaged over settings of the
  % hyperparameters; (2) for each SNP, the posterior inclusion probability
  % averaged over settings of the hyperparameters.
  fprintf('Loading null hypothesis results.\n');
  null            = load(strcat(resultsdir,'/ra/pathway-ra-D.mat'));
  w0              = normalizelogweights(null.logw);
  PIP0            = null.alpha * w0;
  [P1null P2null] = compilesegstats(Ar,null.alpha,w0);

  % ENRICHMENT OF "MEASLES" + MHC
  % Load the results conditioned on the hypothesis that the "Measles"
  % pathway is enriched.
  fprintf('Loading pathway enrichment results.\n');
  load(strcat(resultsdir,'/ra/pathway-ra-J.mat'));
  i         = 3;
  n1        = numel(theta);
  logw1     = logw1(:,:,i);
  alphapath = alpha{i};
  mupath    = mu{i};
  spath     = s{i};

  % Re-load the results conditioned on the null hypothesis that the MHC is
  % enriched, this time selecting the hyperparameter settings that match the
  % hyperparameter settings under the alternative hypothesis.
  null = getnullstats(strcat(resultsdir,'/ra/pathway-ra-D.mat'),theta0);
  mhcnull = load(strcat(resultsdir,'/ra/pathway-ra-G.mat'));
  [logw0 alpha0 mu0 s0] = getpathwaynull(null,mhcnull,theta0,thetaMHC,1);

  % Get the SNPs that lie within the MHC region.
  p = length(labels);
  isinmhc = zeros(p,1);
  isinmhc(mhcsnps) = 1;

  % Get the index of the enriched pathway (j), and the set of SNPs assigned to
  % the enriched pathway that are not in the MHC (snps). Then get the
  % posterior inclusion probabilities and other posterior expectations for
  % each setting of the hyperparameters.
  j     = find(H(:,i));
  snps  = find(A(:,j) & ~isinmhc);
  alpha = repmat(alpha0,1,n1);
  mu    = repmat(mu0,1,n1);
  s     = repmat(s0,1,n1);

  alpha(snps,:) = alphapath(snps,:);
  mu(snps,:)    = mupath(snps,:);
  s(snps,:)     = spath(snps,:);

  % Compute the posterior inclusion probabilities (PIPs) averaged over the
  % settings of the hyperparameters and, for each selected region, compute
  % the posterior statistics P1 and P2 averaged over settings of the
  % hyperparameters.
  w       = normalizelogweights(logw1);
  PIP     = alpha * w(:);
  [P1 P2] = compilesegstats(Ar,alpha,w(:));

  % Show the table legend.
  fprintf('RESULTS (PART O)\n');
  fprintf('                        P1   P1   P2   P2         ');
  fprintf('                                  MAF   MAF\n');
  fprintf('dis chr  region (Mb)  null  alt null  alt gene(s) ');
  fprintf('SNP         PIP   LOR(95%% CI)    ctrls cases\n');

  % Repeat for each selected region.
  n = length(region.label);
  for i = 1:n

    % Get the SNPs in the region.
    snps = find(Ar(:,i));
          
    % Get the SNP in the region with the largest posterior inclusion
    % probability.
    [ans k] = max(PIP(snps));
    k       = snps(k);
      
    % Get the mean and 95% credible interval for the additive effect
    % (i.e. log-odds ratio) corresponding to the SNP with the largest
    % PIP in the segment.
    E     = dot(w(:),mu(k,:));
    beta  = (-2:0.01:2)';
    r     = normpdf(repmat(beta,1,numel(w)),...
                    repmat(mu(k,:),length(beta),1),...
                    repmat(s(k,:),length(beta),1));
    [a b] = cred(beta,r*w(:),E);

    % Print statistics about the selected region.
    fprintf(' RA %3d %6.2f-%-6.2f %0.2f %0.2f %0.2f %0.2f %7s ',...
            region.chr(i),region.start(i)/1e6,region.stop(i)/1e6,...
            P1null(i),P1(i),P2null(i),P2(i),region.label{i});
    fprintf('%-10s %0.2f ',sprintf('rs%d',labels(k)),PIP(k));
    fprintf('%+0.2f(%0.2f-%0.2f) %0.3f %0.3f\n',...
            E,min(abs([a b])),max(abs([a b])),...
            maf(X(ctrls,k)),maf(X(cases,k)));
  end

 case 'P'

  % SCATTERPLOT OF T1D ASSOCIATIONS FOR IL-2 SIGNALING
  % --------------------------------------------------
  % These are labels corresponding to selected segments.
  seglabel = cell(2e4,1);
  seglabel{337}   = 'JUN';
  seglabel{382}   = 'JAK2';
  seglabel{705}   = 'PTPN22';
  seglabel{1133}  = 'MAPKAPK2';
  seglabel{1724}  = 'SOS1';
  seglabel{1781}  = 'PRKCE';
  seglabel{2384}  = 'STAM2';
  seglabel{2610}  = 'STAT1';
  seglabel{2832}  = 'IRS1';
  seglabel{3051}  = 'RAF1';
  seglabel{3286}  = 'CISH';
  seglabel{4044}  = 'PIK3CA';
  seglabel{4902}  = 'IL2';
  seglabel{5742}  = 'PIK3R1';
  seglabel{5850}  = 'RASA1';
  seglabel{6467}  = 'MAPK9';
  seglabel{6761}  = 'MAPK14';
  seglabel{7217}  = 'FYN';
  seglabel{8776}  = 'DOK2';
  seglabel{8818}  = 'PTK2B';
  seglabel{9404}  = 'MYC';
  seglabel{10047} = 'SYK';
  seglabel{10435} = 'IL2RA';
  seglabel{10549} = 'STAM';
  seglabel{10743} = 'MAPK8';
  seglabel{11901} = 'GAB2';
  seglabel{12720} = 'ERBB3(CDK2)';
  seglabel{12797} = 'IFNG';
  seglabel{12957} = 'SOCS2';
  seglabel{13092} = 'SH2B3(PTPN11)';
  seglabel{13910} = 'IRS2';
  seglabel{14299} = 'FOS';
  seglabel{14786} = 'MAP2K1';
  seglabel{15117} = 'CLEC16A(SOCS1)';
  seglabel{15188} = 'PRKCB';
  seglabel{15740} = 'ORMDL3(IKZF3)';
  seglabel{15753} = 'STAT3';
  seglabel{15934} = 'GRB2';
  seglabel{16369} = 'BCL2';
  seglabel{17458} = 'MAPK1';
  seglabel{17572} = 'C1QTNF6(IL2RB)';

  figure(1);
  clf
  windowpos = get(gcf,'Position');
  set(gcf,'Name','Results (Part P)','Color','white',...
          'NumberTitle','off','PaperPositionMode','auto',...
          'Position',[windowpos(1:2) 1050 425]);

  % Load the genotype and phenotype data.
  fprintf('Loading genotype and phenotype data for T1D.\n');
  load(strcat(datadir,'/t1d.mat'));

  % Create overlapping segments with 50 SNPs.
  fprintf('Creating segments.\n');
  [Aseg segchr segpos] = makesegments2(chr,pos,25,5e4);
  [Aseg segchr segpos] = combinesegments(Aseg,segchr,segpos);

  % NULL HYPOTHESIS - NO ENRICHED PATHWAYS
  % Load the results conditioned on the null hypothesis that no pathways are
  % enriched and, for each segment, compute the posterior statistics P1 and
  % P2 averaged over settings of the hyperparameters.
  fprintf('Loading null hypothesis results.\n');
  null            = load(strcat(resultsdir,'/t1d/pathway-t1d-D.mat'));
  w0              = normalizelogweights(null.logw);
  [P1null P2null] = compilesegstats(Aseg,null.alpha,w0);
  
  % ENRICHMENT OF IL-2 SIGNALING + MHC
  % Load the results conditioned on the hypothesis that the IL-2 signaling
  % pathway and SNPs within the MHC region are enriched.
  fprintf('Loading pathway enrichment results.\n');
  load(strcat(resultsdir,'/t1d/pathway-t1d-O.mat'));
  i         = 8;
  n1        = numel(theta);
  logw1     = logw1(:,:,i);
  alphapath = alpha{i};
  mupath    = mu{i};

  % Re-load the results conditioned on the null hypothesis that the MHC is
  % enriched, this time selecting the hyperparameter settings that match the
  % hyperparameter settings under the alternative hypothesis.
  null = getnullstats(strcat(resultsdir,'/t1d/pathway-t1d-D.mat'),theta0);
  mhcnull = load(strcat(resultsdir,'/t1d/pathway-t1d-G.mat'));
  [logw0 alpha0 mu0 s0] = getpathwaynull(null,mhcnull,theta0,thetaMHC,1);

  % Get the SNPs that lie within the MHC region.
  p = length(labels);
  isinmhc = zeros(p,1);
  isinmhc(mhcsnps) = 1;

  % Get the index of the enriched pathway (j), and the set of SNPs assigned to
  % the enriched pathway that are not in the MHC (snps). Then get the
  % posterior inclusion probabilities (and posterior mean additive effects)
  % for each setting of the hyperparameters.
  j             = find(H(:,i));
  snps          = find(A(:,j) & ~isinmhc);
  alpha         = repmat(alpha0,1,n1);
  mu            = repmat(mu0,1,n1);
  alpha(snps,:) = alphapath(snps,:);
  mu(snps,:)    = mupath(snps,:);

  % Compute the posterior inclusion probabilities (PIPs) and posterior mean
  % regression coefficients averaged over the settings of the
  % hyperparameters and, for each segment, compute the posterior statistics
  % P1 and P2 averaged over settings of the hyperparameters.
  w       = normalizelogweights(logw1);
  PIP     = alpha * w(:);
  mu      = mu * w(:);
  [P1 P2] = compilesegstats(Aseg,alpha,w(:));

  % Plot the posterior probabilities given that no pathways are enriched
  % (horizontal axis) versus the posterior probabilities given that the IL-2
  % signaling pathway, including the MHC, are enriched (vertical axis). I
  % highlight segments that overlap SNPs assigned to the enriched pathway.
  subplot(1,2,1);
  plot([0 2],[0 2],':','Color',rgb('dodgerblue'),'LineWidth',1);
  hold on
  inpathseg = full(sumrows(Aseg(snps,:)))' > 0;
  inmhcseg  = full(sumrows(Aseg(mhcsnps,:)))' > 0;

  % First, plot "non-redundant" segments overlapping the MHC.
  I = (P1 > 0.01 | P1null > 0.01) & inmhcseg & ~redundantsegment(Aseg,P1);
  plot(P1null(I),P1(I),'o','MarkerFaceColor',rgb('royalblue'),...
       'MarkerEdgeColor','none','MarkerSize',6);

  % Next, plot "non-redundant" segments that have posterior probabilities
  % less than 0.4 under the alternative hypothesis.
  I = (P1 > 0.01 | P1null > 0.01) & P1 < 0.4 & ~inmhcseg& ...
      ~redundantsegment(Aseg,P1);
  t = I & ~inpathseg;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('black'),...
       'MarkerEdgeColor','none','MarkerSize',6);
  t = I & inpathseg;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('firebrick'),...
       'MarkerEdgeColor','none','MarkerSize',6);

  % Finally, plot "non-redundant" segments with posterior probabilities
  % greater than 0.4 under the alternative hypothesis. Here I need to fix a
  % few segments that underestimate the posterior probability under the null
  % hypothesis.
  I            = P1 >= 0.4 & ~inmhcseg & ~redundantsegment(Aseg,P1);
  I(15117)     = 1;
  I(15118)     = 0;
  P1null(5850) = P1null(5849);
  t = I & ~inpathseg;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('black'),...
       'MarkerEdgeColor','none','MarkerSize',6);
  t = I & inpathseg;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('firebrick'),...
       'MarkerEdgeColor','none','MarkerSize',6);

  % Label some of the segments.
  t = I & ~inpathseg;
  text(P1null(t) + 0.005,P1(t) + 0.005,seglabel(t,:),...
       'VerticalAlignment','bottom','FontSize',9,'FontName','fixed',...
       'Color',rgb('black'));
  t = I & inpathseg;
  text(P1null(t) + 0.005,P1(t) + 0.005,seglabel(t,:),...
       'VerticalAlignment','bottom','FontSize',9,'FontName','fixed',...
       'Color',rgb('firebrick'));

  hold off
  set(gca,'FontSize',10,'FontName','fixed');
  set(gca,'XLim',[0 1.05],'XTick',0:0.5:1);
  set(gca,'YLim',[0 1.05],'YTick',0:0.5:1);
  set(gca,'TickDir','out');
  xlabel('P1 when no pathways are enriched');
  ylabel('P1 when MHC and IL-2 signaling pathway are enriched');
  title('Type 1 diabetes');

  % SCATTERPLOT OF T1D ASSOCIATIONS FOR 2 ENRICHED PATHWAYS
  % -------------------------------------------------------
  % Save the results for enrichment of IL-2 signaling + MHC.
  P1null = P1;
  P2null = P2;

  % Load the results conditioned on enrichment of 2 pathways.
  fprintf('ENRICHMENT OF 2 PATHWAYS.\n');
  fprintf('Loading pathway enrichment results.\n');
  load(strcat(resultsdir,'/t1d/pathway-t1d-L.mat'));
  n1        = numel(theta);
  alphapath = alpha;

  % Re-load the results conditioned on the null hypothesis that the MHC is
  % enriched, this time selecting the hyperparameter settings that match the
  % hyperparameter settings under the alternative hypothesis.
  null = getnullstats(strcat(resultsdir,'/t1d/pathway-t1d-D.mat'),theta0);
  mhcnull = load(strcat(resultsdir,'/t1d/pathway-t1d-G.mat'));
  [logw0 alpha0 mu0 s0] = getpathwaynull(null,mhcnull,theta0,thetaMHC,1);

  % Compute the posterior statistics P1 for each enrichment hypothesis.
  n    = length(BF);
  nseg = length(segchr);
  P1   = zeros(nseg,n);
  P2   = zeros(nseg,n);
  fprintf('Computing segment statistics given enrichment hypotheses:\n');
  for i = 1:n

    % Get the indices of the enriched pathways (paths), and the set of SNPs
    % assigned to the enriched pathways (snps).
    paths = find(H(:,i))';
    snps  = find(sumcols(A(:,paths)));
    fprintf('%d. Enrichment of MHC + "%s" + "%s"\n',i,...
            pathway.label{paths(1)},...
            pathway.label{paths(2)});

    % Get the posterior inclusion probabilities for each setting of the
    % hyperparameters.
    alpha         = repmat(alpha0,[1 n1]);
    alpha(snps,:) = alphapath{i}(snps,:);

    % For each segment, compute the posterior statistics P1 and P2
    % averaged over settings of the hyperparameters.
    w                 = normalizelogweights(logw1(:,:,i));
    [P1(:,i) P2(:,i)] = compilesegstats(Aseg,alpha,w(:));      
  end

  % Get the final posterior statistics P1 and P2 by averaging over the
  % enrichment hypotheses.
  P1 = P1 * BF / sum(BF);
  P2 = P2 * BF / sum(BF);

  % Plot the posterior probabilities given that the IL-2 signaling and the MHC
  % region are enriched (horizontal axis) versus the posterior probabilities
  % given that 2 pathways plus the MHC are enriched (vertical axis). I
  % highlight segments that overlap SNPs assigned to the additional enriched
  % pathways (and are not assigned to the IL-2 signaling pathway, nor to the
  % MHC).
  j         = 2608;
  paths     = find(sumcols(H));
  snps      = find(sumcols(A(:,paths)) & ~A(:,j) & ~isinmhc);
  inpathseg = full(sumrows(Aseg(snps,:)))' > 0;
  inmhcseg  = full(sumrows(Aseg(mhcsnps,:)))' > 0;
  subplot(1,2,2);
  plot([0 2],[0 2],':','Color',rgb('dodgerblue'),'LineWidth',1);
  hold on

  % First, plot "non-redundant" segments overlapping the MHC.
  I = (P1 > 0.01 | P1null > 0.01) & inmhcseg & ~redundantsegment(Aseg,P1);
  plot(P1null(I),P1(I),'o','MarkerFaceColor',rgb('royalblue'),...
       'MarkerEdgeColor','none','MarkerSize',6);

  % Plot "non-redundant" segments that do not overlap the MHC.
  I = (P1 > 0.01 | P1null > 0.01) & ~inmhcseg & ~redundantsegment(Aseg,P1);
  t = I & ~inpathseg;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('black'),...
       'MarkerEdgeColor','none','MarkerSize',6);
  t = I & inpathseg;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('firebrick'),...
       'MarkerEdgeColor','none','MarkerSize',6);

  % Label some of the segments.
  t = I & inpathseg;
  text(P1null(t) + 0.005,P1(t) + 0.005,seglabel(t,:),...
       'VerticalAlignment','bottom','FontSize',9,...
       'FontName','fixed','Color',rgb('firebrick'));

  hold off
  set(gca,'FontSize',10,'FontName','fixed');
  set(gca,'XLim',[0 1.05],'XTick',0:0.1:1);
  set(gca,'YLim',[0 1.05],'YTick',0:0.1:1);
  set(gca,'TickDir','out');
  xlabel('P1 when IL-2 signaling + MHC are enriched');
  ylabel('P1 when 2 pathways + MHC are enriched');
  title('type 1 diabetes');

 case 'Q'

  % SELECTED REGIONS OF GENOME FOR T1D GIVEN ENRICHMENT OF IL-2 SIGNALING
  % ---------------------------------------------------------------------
  % The columns of this table specify: (1-3) the chromosome and the start and
  % end positions on the chromosome (in Mb) that make up the "critical
  % region" for the disease association; (4) candidate gene(s) corresponding
  % to locus.
  D = {  1 113.53 114.36 'PTPN22' 
         3 12.44   12.85 'RAF1'
         4 123.30 123.93 'IL2'
         6  25.52  33.76 'MHC'
         6  35.86  36.38 'MAPK14'
         6 112.02 112.47 'FYN'
        10   6.07   6.26 'IL2RA'
        12  54.63  54.91 'ERBB3'
        12 109.80 111.74 'SH2B3' 
        16  10.90  11.41 'CLEC16A'
        17  34.60  35.50 'ORMDL3'
        22  35.68  36.00 'C1QTNF6' };

  region.chr   = cell2mat(D(:,1));
  region.start = cell2mat(D(:,2)) * 1e6;
  region.stop  = cell2mat(D(:,3)) * 1e6;
  region.label = D(:,4);

  % Load the genotype and phenotype data.
  fprintf('Loading genotype and phenotype data for T1D.\n');
  load(strcat(datadir,'/t1d.mat'));

  % Get the cases and controls.
  cases = find(y == 1);
  ctrls = find(y == 0);

  % Get the set of SNPs for each selected region.
  Ar = snps2genes(chr,pos,region,0,1e4);

  % NULL HYPOTHESIS - NO ENRICHED PATHWAYS
  % Load the results conditioned on the null hypothesis that no pathways are
  % enriched and compute: (1) for each of the selected regions, the
  % posterior statistics P1 and P2 averaged over settings of the
  % hyperparameters; (2) for each SNP, the posterior inclusion probability
  % averaged over settings of the hyperparameters.
  fprintf('Loading null hypothesis results.\n');
  null            = load(strcat(resultsdir,'/t1d/pathway-t1d-D.mat'));
  w0              = normalizelogweights(null.logw);
  PIP0            = null.alpha * w0;
  [P1null P2null] = compilesegstats(Ar,null.alpha,w0);

  % ENRICHMENT OF IL-2 SIGNALING + MHC
  % Load the results conditioned on the hypothesis that the IL-2 signaling
  % pathway is enriched.
  fprintf('Loading pathway enrichment results.\n');
  load(strcat(resultsdir,'/t1d/pathway-t1d-O.mat'));
  i         = 8;
  n1        = numel(theta);
  logw1     = logw1(:,:,i);
  alphapath = alpha{i};
  mupath    = mu{i};
  spath     = s{i};

  % Re-load the results conditioned on the null hypothesis that the MHC is
  % enriched, this time selecting the hyperparameter settings that match the
  % hyperparameter settings under the alternative hypothesis.
  null = getnullstats(strcat(resultsdir,'/t1d/pathway-t1d-D.mat'),theta0);
  mhcnull = load(strcat(resultsdir,'/t1d/pathway-t1d-G.mat'));
  [logw0 alpha0 mu0 s0] = getpathwaynull(null,mhcnull,theta0,thetaMHC,1);

  % Get the SNPs that lie within the MHC region.
  p = length(labels);
  isinmhc = zeros(p,1);
  isinmhc(mhcsnps) = 1;

  % Get the index of the enriched pathway (j), and the set of SNPs assigned to
  % the enriched pathway that are not in the MHC (snps). Then get the
  % posterior inclusion probabilities and other posterior expectations for
  % each setting of the hyperparameters.
  j     = find(H(:,i));
  snps  = find(A(:,j) & ~isinmhc);
  alpha = repmat(alpha0,1,n1);
  mu    = repmat(mu0,1,n1);
  s     = repmat(s0,1,n1);

  alpha(snps,:) = alphapath(snps,:);
  mu(snps,:)    = mupath(snps,:);
  s(snps,:)     = spath(snps,:);
  
  % Compute the posterior inclusion probabilities (PIPs) averaged over the
  % settings of the hyperparameters and, for each selected region, compute
  % the posterior statistics P1 and P2 averaged over settings of the
  % hyperparameters.
  w       = normalizelogweights(logw1);
  PIP     = alpha * w(:);
  [P1 P2] = compilesegstats(Ar,alpha,w(:));
  
  % Show the table legend.
  fprintf('RESULTS (PART Q)\n');
  fprintf('                        P1   P1   P2   P2         ');
  fprintf('                                  MAF   MAF\n');
  fprintf('dis chr  region (Mb)  null  alt null  alt gene(s) ');
  fprintf('SNP         PIP   LOR(95%% CI)    ctrls cases\n');

  % Repeat for each selected region.
  n = length(region.label);
  for i = 1:n

    % Get the SNPs in the region.
    snps = find(Ar(:,i));
          
    % Get the SNP in the region with the largest posterior inclusion
    % probability.
    [ans k] = max(PIP(snps));
    k       = snps(k);
      
    % Get the mean and 95% credible interval for the additive effect
    % (i.e. log-odds ratio) corresponding to the SNP with the largest
    % PIP in the segment.
    E     = dot(w(:),mu(k,:));
    beta  = (-2:0.01:2)';
    r     = normpdf(repmat(beta,1,numel(w)),...
                    repmat(mu(k,:),length(beta),1),...
                    repmat(s(k,:),length(beta),1));
    [a b] = cred(beta,r*w(:),E);

    % Print statistics about the selected region.
    fprintf('T1D %3d %6.2f-%-6.2f %0.2f %0.2f %0.2f %0.2f %7s ',...
            region.chr(i),region.start(i)/1e6,region.stop(i)/1e6,...
            P1null(i),P1(i),P2null(i),P2(i),region.label{i});
    fprintf('%-10s %0.2f ',sprintf('rs%d',labels(k)),PIP(k));
    fprintf('%+0.2f(%0.2f-%0.2f) %0.3f %0.3f\n',...
            E,min(abs([a b])),max(abs([a b])),...
            maf(X(ctrls,k)),maf(X(cases,k)));
  end

 case 'R'

  % SCATTERPLOT OF T2D ASSOCIATIONS WITH/WITHOUT PATHWAYS
  % -----------------------------------------------------
  % These are labels corresponding to the selected segments.
  seglabel = cell(2e4,1);
  seglabel{11077} = 'GPR120';
  seglabel{11220} = 'TCF7L2';
  seglabel{15260} = 'FTO';
  
  figure(1);
  clf
  windowpos = get(gcf,'Position');
  set(gcf,'Name','Results (Part R)','Color','white',...
          'NumberTitle','off','PaperPositionMode','auto',...
          'Position',[windowpos(1:2) 640 570]);

  % Load the genotype and phenotype data.
  fprintf('Loading genotype and phenotype data for T2D.\n');
  load(strcat(datadir,'/t2d.mat'));

  % Create overlapping segments with 50 SNPs.
  fprintf('Creating segments.\n');
  [Aseg segchr segpos] = makesegments2(chr,pos,25,5e4);
  [Aseg segchr segpos] = combinesegments(Aseg,segchr,segpos);

  % NULL HYPOTHESIS - NO ENRICHED PATHWAYS
  % Load the results conditioned on the null hypothesis that no pathways are
  % enriched and, for each segment, compute the posterior statistics P1 and
  % P2 averaged over settings of the hyperparameters.
  fprintf('Loading null hypothesis results.\n');
  null            = load(strcat(resultsdir,'/t2d/pathway-t2d-B.mat'));
  w0              = normalizelogweights(null.logw);
  [P1null P2null] = compilesegstats(Aseg,null.alpha,w0);

  % ENRICHMENT HYPOTHESES WITH 1 ENRICHED PATHWAY
  % Load the results conditioned on hypotheses in which 1 pathway is enriched.
  fprintf('Loading pathway enrichment results.\n');
  load(strcat(resultsdir,'/t2d/pathway-t2d-E.mat'));
  n1        = numel(theta);
  alphapath = alpha;

  % Re-load the results conditioned on the null hypothesis that no pathways
  % are enriched, this time selecting the hyperparameter settings that match
  % the hyperparameter settings under the alternative hypothesis.
  null = getnullstats(strcat(resultsdir,'/t2d/pathway-t2d-B.mat'),theta0);

  % Compute the posterior statistics P1 and P2 for each enrichment hypothesis.
  n    = length(BF);
  nseg = length(segchr);
  P1   = zeros(nseg,n);
  P2   = zeros(nseg,n);
  fprintf('Computing segment statistics given enrichment hypotheses:\n');
  for i = 1:n
    
    % Get the index of the enriched pathway (j), and the set of SNPs
    % assigned to the enriched pathway (snps).
    j    = find(H(:,i));
    snps = find(A(:,j));
    fprintf('%d. Enrichment of "%s"\n',i,pathway.label{j});

    % Get the posterior inclusion probabilities for each setting of the
    % hyperparameters.
    alpha         = repmat(null.alpha,1,n1);
    alpha(snps,:) = alphapath{i}(snps,:);

    % For each segment, compute the posterior statistics P1 and P2
    % averaged over settings of the hyperparameters.
    w                 = normalizelogweights(logw1(:,:,i));
    [P1(:,i) P2(:,i)] = compilesegstats(Aseg,alpha,w(:));
  end

  % Get the final posterior statistics P1 and P2 by averaging over the
  % enrichment hypotheses.
  P1 = P1 * BF / sum(BF);
  P2 = P2 * BF / sum(BF);
  
  % Plot the posterior probabilities P1 given that no pathways are enriched
  % (horizontal axis) versus the posterior probabilities P1 given that one
  % pathway is enriched (vertical axis). I highlight segments that overlap
  % SNPs assigned to the enriched pathway. There is a segment near genes
  % RBL2 and FTO that is assigned to one of the enriched pathways---Id
  % signaling pathway---but does not appear to have much effect on the
  % posterior probability P1, so I do not show it as overlapping with SNPs
  % in the pathway.
  paths  = find(sumcols(H));
  snps   = find(sumcols(A(:,paths)));
  inpath = full(sumrows(Aseg(snps,:)))' > 0;
  inpath(15260) = 0;
  plot([0 2],[0 2],':','Color',rgb('dodgerblue'),'LineWidth',1);
  hold on

  % Plot "non-redundant" segments.
  I = (P1 > 0.01 | P1null > 0.01) & ~redundantsegment(Aseg,P1);
  t = I & ~inpath;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('black'),...
       'MarkerEdgeColor','none','MarkerSize',6);
  t = I & inpath;
  plot(P1null(t),P1(t),'o','MarkerFaceColor',rgb('firebrick'),...
       'MarkerEdgeColor','none','MarkerSize',6);

  % Label some of the segments.
  t = I & ~inpath;
  text(P1null(t) + 0.005,P1(t) + 0.005,seglabel(t,:),...
       'VerticalAlignment','bottom','FontSize',9,...
       'FontName','fixed','Color',rgb('black'));
  t = I & inpath;
  text(P1null(t) + 0.005,P1(t) + 0.005,seglabel(t,:),...
       'VerticalAlignment','bottom','FontSize',9,...
       'FontName','fixed','Color',rgb('firebrick'));

  hold off
  set(gca,'FontSize',10,'FontName','fixed');
  set(gca,'XLim',[0 1.05],'XTick',0:0.1:1);
  set(gca,'YLim',[0 1.05],'YTick',0:0.1:1);
  set(gca,'TickDir','out');
  xlabel('P1 when no pathways are enriched');
  ylabel('P1 when 1 pathway is enriched');
  title('type 2 diabetes');

 case 'S'

  % BAYES FACTORS FOR COMBINATIONS OF ENRICHED PATHWAYS
  % ---------------------------------------------------
  % The columns of this table specify: (1) the disease; (2) the file
  % containing the more accurate estimates of the Bayes factors; (3) if this
  % number is greater than 1, then the Bayes factor is calculated with
  % respect to the hypothesis that the MHC is enriched for disease risk
  % factors, in which case we need to multiply by this number to obtain the
  % Bayes factor with respect to the hypothesis that no pathways are
  % enriched; (4) the combination of pathways included in the top enrichment
  % hypothesis, and the order in which to display these pathways in the
  % results.
  D = { 't1d' 't1d/pathway-t1d-N.mat' 2.70e54 [3322 2608 2576 1136]
        't1d' 't1d/pathway-t1d-L.mat' 2.70e54 [3322 2608 2576]
        't1d' 't1d/pathway-t1d-J.mat' 2.70e54 [3322 2608]
        't1d' 't1d/pathway-t1d-G.mat' 1       3322
        'ra'  'ra/pathway-ra-N.mat'   2.43e21 [3322 2088 1447 1136]
        'ra'  'ra/pathway-ra-L.mat'   2.43e21 [3322 2088 1447]
        'ra'  'ra/pathway-ra-J.mat'   2.43e21 [3322 2088]
        'ra'  'ra/pathway-ra-G.mat'   1       3322
        'cd'  'cd/pathway-cd-J.mat'   1       [2069 2591 684]
        'cd'  'cd/pathway-cd-G.mat'   1       [2069 2591]
        'cd'  'cd/pathway-cd-D.mat'   1       2069 };

  disease       = D(:,1);
  resultsfile   = D(:,2);
  BFmhc         = cell2mat(D(:,3));
  pathways      = D(:,4);

  figure(1);
  clf
  windowpos = get(gcf,'Position');
  set(gcf,'Name','Results (Part S)','Color','white',...
          'NumberTitle','off','PaperPositionMode','auto',...
          'Position',[windowpos(1:2) 1300 325]);

  % Initialize storage for the combined results.
  n        = length(disease);
  numgenes = zeros(n,4);
  numsnps  = zeros(n,4);
  E0       = zeros(n,1);
  a0       = zeros(n,1);
  b0       = zeros(n,1);
  E        = zeros(n,1);
  a        = zeros(n,1);
  b        = zeros(n,1);

  % Repeat for each set of results.
  for i = 1:n

    % Load the Bayes factors and other statistics.
    load(strcat(resultsdir,'/',resultsfile{i}));

    % Get the total number of genes and SNPs.
    ng = length(gene.symbol);
    p  = size(A,1);

    % Get results for the enrichment hypothesis with the largest Bayes
    % factor. 
    [ans k] = max(BF);
    BF      = BF(k);
    H       = H(:,k);
    logw1   = logw1(:,:,k);

    % Get the number of genes (and SNPs) assigned to each pathway.
    m     = length(pathways{i});
    asnp  = zeros(p,1);
    agene = zeros(ng,1);
    for t = 1:m
      j = pathways{i}(t);

      % Get the number of genes assigned to the pathway that aren't
      % already assigned to other pathways.
      I             = pathway.genes(:,j) & inref(gene) & ~agene;
      numgenes(i,t) = sum(I);
      agene(I)      = 1;

      % Get the number of SNPs assigned to the pathway that aren't
      % already assigned to other pathways.
      I            = A(:,j) & ~asnp;
      numsnps(i,t) = sum(I);
      asnp(I)      = 1;
    end

    % Get the posterior mean and 95% credible interval interval for the
    % genome-wide log-odds (theta0).
    w             = normalizelogweights(logw1);
    E0(i)         = dot(theta0,sumcols(w));
    [a0(i) b0(i)] = cred(theta0,sumcols(w),E0(i));

    % Get the posterior mean and 95% credible interval interval for the
    % enrichment parameter (theta)
    E(i)        = dot(theta,sumrows(w));
    [a(i) b(i)] = cred(theta,sumrows(w),E(i));

    % Print more details about the enrichment hypothesis.
    if m == 1
      fprintf('%s -- 1 enriched pathway:\n',upper(disease{i}));
    else
      fprintf('%s -- %d enriched pathways:\n',upper(disease{i}),m);
    end
    for t = 1:m
      j = pathways{i}(t);
      fprintf('%d. %s (%s/%s,i=%d)\n',t,pathway.label{j},...
              pathway.source{j},pathway.database{j},j);
    end
    BF = BF * BFmhc(i);
    if BF > 1e4
      fprintf('BF     = %0.2e\n',BF);
    else
      fprintf('BF     = %d\n',round(BF));
    end
    fprintf('theta0 = %0.2f (%0.1f,%0.1f)\n',a0(i),E0(i),b0(i));
    fprintf('theta  = %0.2f (%0.1f,%0.1f)\n',a(i),E(i),b(i));
    fprintf('\n');
  end

  % Plot the number of genes.
  subplot(1,6,2:3);
  h = barh(1:n,numgenes,0.6,'stacked');
  set(h,'EdgeColor',rgb('white'),'LineWidth',3,...
        'FaceColor',rgb('darkorange'));
  set(gca,'FontSize',9,'FontName','fixed');
  set(gca,'YDir','reverse','YLim',[0 n+1],'YTick',1:n);
  set(gca,'XLim',[0 400],'XTick',0:100:400);
  set(gca,'TickLength',[0 0]);
  xlabel('no. genes','FontSize',10);
  str = cell(n,1);
  for i = 1:n
    m      = length(pathways{i});
    str{i} = num2str(numgenes(i,1));
    if m > 1
      str{i} = strcat(str{i},sprintf('+%d',numgenes(i,2:m)));
    end
  end
  text(sumcols(numgenes) + 5,1:n,str,'FontSize',9,'FontName','fixed');

  % Label the bar plot with the disease and the number of enriched pathways.
  str = cell(n,1);
  for i = 1:n
    m = length(pathways{i});
    if m == 1
      str{i} = sprintf('%s-1 enriched pathway',upper(disease{i}));
    else
      str{i} = sprintf('%s-%d enriched pathways',upper(disease{i}),m);
    end
  end
  set(gca,'YTickLabel',str);

  % Plot the number of SNPs.
  subplot(1,6,4);
  h = barh(1:n,numsnps,0.6,'stacked');
  set(h,'EdgeColor',rgb('white'),'LineWidth',3,...
        'FaceColor',rgb('darkorange'));
  set(gca,'FontSize',9,'FontName','fixed');
  set(gca,'YDir','reverse','YLim',[0 n+1],'YTick',[]);
  set(gca,'XLim',[0 1e4],'XTick',[0 1e4]);
  xlabel('no. SNPs','FontSize',10);
  str = cell(n,1);
  for i = 1:n
    m      = length(pathways{i});
    str{i} = num2str(numsnps(i,1));
    if m > 1
      str{i} = strcat(str{i},sprintf('+%d',numsnps(i,2:m)));
    end
  end
  text(sumcols(numsnps) + 200,1:n,str,'FontSize',9,'FontName','fixed');

  % Show the posterior mean and 95% credible interval of the genome-wide
  % log-odds (theta0) for each enrichment hypothesis.
  subplot(1,6,5);
  h = herrorbar(E0,1:n,E0-a0,b0-E0,'ko');
  set(h,'Markersize',4,'MarkerFaceColor',rgb('black'));
  set(gca,'FontSize',9,'FontName','fixed');
  set(gca,'YDir','reverse','YLim',[0 n+1],'YTick',[]);
  set(gca,'XLim',[-6.25 -2],'XTick',-6:-2,'XGrid','on');
  set(gca,'TickLength',[0 0]);
  xlabel('genome-wide log-odds','FontSize',10);

  % Show the posterior mean and 95% credible interval of the enrichment
  % parameter (theta) for each enrichment hypothesis.
  subplot(1,6,6);
  h = herrorbar(E,1:n,E-a,b-E,'ko');
  set(h,'Markersize',4,'MarkerFaceColor',rgb('black'));
  set(gca,'FontSize',9,'FontName','fixed');
  set(gca,'YDir','reverse','YLim',[0 n+1],'YTick',[]);
  set(gca,'XLim',[0 5],'XTick',0:5,'XGrid','on');
  set(gca,'TickLength',[0 0]);
  xlabel('log_{10}enrichment','FontSize',10);
end
