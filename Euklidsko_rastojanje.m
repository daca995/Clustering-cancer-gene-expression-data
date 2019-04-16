function [E]= Euklidsko_rastojanje(X)

%X = importdata('armstrong_v1.txt');

for i = 1:size(X,1)
    for j = 1:size(X,1)
        E(i,j) = sqrt(sum((X(i,:) - X(j,:)).^2));
    end
end
