% INMHCREGION(CHR,POS) returns a logical vector of the same size as inputs
% CHR and POS that is equal to TRUE if the SNP is within the "extended"
% major histocompatibility complex (MHC). Input CHR specifies the chromosome
% number of each SNP, and POS gives the position on the chromosome.
function I = inmhcregion (chr, pos)
  I = chr == 6 & pos > 25.7e6 & pos < 33.6e6;

                     