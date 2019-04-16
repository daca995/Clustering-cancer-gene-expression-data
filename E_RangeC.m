function [R] = E_RangeC(M)
     M = (M- min(M))./(max(M)-min(M));
     R = Euklidsko_rastojanje(M);
end
