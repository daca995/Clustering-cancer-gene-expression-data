function [R] = E_ZscoreC(M)
     M = zscore(M);
     R = Euklidsko_rastojanje(M);
end
