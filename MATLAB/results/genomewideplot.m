% GENOMEWIDEPLOT(CHR,POS,SCORE) shows the genome-wide scan (the "Manhattan
% plot" in the WTCCC paper). Inputs CHR and POS specify the chromosome
% numbers and chromosomal positions of the markers. The SNPs are assumed to
% be ordered by chromosome number, then by position along the chromosome.
% SCORE is some measure of significance for all the markers, such as Bayes
% factor, posterior odds, or posterior inclusion probability.
function genomewideplot (chr, pos, score, snps)

  % Alternating shades for chromosomes.
  clr1 = rgb('midnightblue');
  clr2 = rgb('cornflowerblue');

  % These are the SNPs to highlight in the Manhattan plot.
  if ~exist('snps')
    snps = [];
  end
  
  % Get the set of chromosomes represented by the SNPs.
  chrs = unique(chr);
  chrs = chrs(:)';
  
  % This is the current position in the plot.
  p0 = 0;

  % Repeat for each chromosome.
  hold on
  for c = chrs
    
    % Plot the SNPs on the chromosome.
    is = find(chr == c);
    maxpos = max(pos(is));
    if isodd(c)
      clr = clr1;
    else
      clr = clr2;
    end
    plot(p0 + pos(is),score(is),'o','MarkerFaceColor',clr,...
	 'MarkerEdgeColor','none','MarkerSize',4);

    % Highlight these SNPs.
    is = find(chr(snps) == c);
    if isodd(c)
      clr = rgb('forestgreen');
    else
      clr = rgb('yellowgreen');
    end
    if length(is)
      plot(p0 + pos(snps(is)),score(snps(is)),'o','MarkerFaceColor',clr,...
	   'MarkerEdgeColor','none','MarkerSize',4);    
    end

    % Add the chromosome number.
    text(p0 + maxpos/2,-0.25 * (max(score) - min(score)),num2str(c),...
	 'Color','black','FontSize',9,'HorizontalAlignment','center',...
	 'VerticalAlignment','bottom','FontName','Arial');

    % Move to the next chromosome.
    p0 = p0 + maxpos;
  end
  hold off
  set(gca,'XLim',[0 p0],'XTick',[],'YLim',[min(score) max(score)]);
  set(gca,'TickDir','out','TickLength',[0.03 0.03]);

  