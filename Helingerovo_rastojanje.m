function [H]= Helingerovo_rastojanje(X)

%X = importdata('armstrong_v1.txt');
X = sqrt(X./(sum(X'))');

for i = 1:size(X,1)
    for j = 1:size(X,1)
        H(i,j) = sqrt(sum((X(i,:) - X(j,:)).^2))./sqrt(2);
    end
end



