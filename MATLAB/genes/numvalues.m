% NUMVALUES(MAP) returns the number of unique values of a Map object (see
% help for CONTAINERS.MAP for more information).
function n = numvalues (map)
  n = length(unique(cell2mat(map.values)));
  