function [R] = E_RangeR(M)
     M = M';
     M = (M- min(M))./(max(M)-min(M));
     R = Euklidsko_rastojanje(M');
     
end
