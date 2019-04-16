function [R] = E_ZscoreR(M)
     M = zscore(M');
     R = Euklidsko_rastojanje(M');
end
