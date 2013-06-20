% Script to compute the genetic distances between SNPs. Note that the
% SNPs in the HapMap data and those that were typed in the study must be
% sorted by chromosome number, then by position along the chromosome.
clear

% Load the SNP data.
load('cd.mat');
p = length(labels);

% Load the genetic distances for the HapMap genetic markers.
load('hapmap.mat');
hapmap.chr  = [ hapmap.chr;  0   ];
hapmap.pos  = [ hapmap.pos;  Inf ];
hapmap.dist = [ hapmap.dist; 0   ];
hapmap.freq = [ hapmap.freq; 0   ];

% Initialize storage for the genetic distances (in centimorgans, or cM).
dist = zeros(p,1);

% Repeat for each chromosome.
i = 1;  % HapMap locus.
k = 1;  % Study locus.
for c = 1:22
  fprintf('-Chromosome %d\n',c);
  
  % Move to the first HapMap SNP on the chromosome.
  while hapmap.chr(i) < c
    i = i + 1;
  end
  
  % For any SNPs that are located before the first HapMap marker, we have
  % no recombination rate info, so set the genetic distance to zero.
  while pos(k) < hapmap.pos(i)
    dist(k) = 0;
    k = k + 1;
  end
  
  % Repeat for each study SNP on the chromosome.
  while chr(k) == c

    % If the chromosomal position of the current study SNP is greater than
    % the chromosomal position of the next HapMap marker on the same
    % chromosome, move to the next HapMap marker.
    if pos(k) >= hapmap.pos(i+1) & hapmap.chr(i+1) == c
      i = i + 1;
    else
    
      % Calculate the genetic distance from the first HapMap marker to
      % the current SNP.
      rate    = hapmap.freq(i) / 1e6;
      dist(k) = hapmap.dist(i) + rate * (pos(k) - hapmap.pos(i));
      
      % Move to the next SNP. Stop if we've reached the last SNP.
      k = k + 1;
      if k > p
	break
      end
    end
  end
end

% Save the genetic distances.
save('cdgendist.mat','dist','-v7.3');
