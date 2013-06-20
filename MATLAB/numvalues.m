% NUMVALUES(MAP) returns the number of unique values 
function n = numvalues (map)
  n = length(unique(cell2mat(map.values)));
  