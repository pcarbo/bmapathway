% REDUNDANTSEGMENT(A,R) returns a binary vector, in which an entry of the
% vector is 1 if and only if the segment is "redundant." By this, we mean
% that the segment overlaps with another segment that is more highly
% ranked. Input A is a P x N matrix representing assignments of SNPs to
% segments, where P is the number of SNPs, and N is the number of
% segments. Input R is a vector of length R containing the ranking of the
% segments (for example, sums of posterior inclusion probabilities).
function y = redundantsegment (A, r)

  % Get the number of segments.
  n = length(r);

  % Compute the adjacency matrix for segments.
  B = spones(A' * A);

  % Start by declaring none of the segments as redundant.
  y = zeros(n,1);

  % Repeat for each segment.
  for i = 1:n
    
    % If the segment shares the association signal with another segment,
    % and the association signal is at least as strong for the other
    % segment, remove this segment.
    if sum(B(:,i) & ~y & r >= r(i)) > 1
      y(i) = 1;
    end
  end
