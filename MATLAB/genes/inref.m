% INREF(GENE) returns a vector of length equal to the number of genes. An
% entry of the vector is TRUE if and only if the gene maps to the reference
% sequence. See function GETREFSEQDATA for more information.
function y = inref (gene)
  y = ~(gene.start == -1 | gene.stop == -1);