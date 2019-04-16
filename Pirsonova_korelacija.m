function [P]= Pirsonova_korelacija(X)

P = corr(X', 'type', 'Pearson');

end

%X = importdata('armstrong_v1.txt');
% 
% for i = 1:size(X,1)
%     for j = 1:size(X,1)
%         x = X(i,:)-mean(X(i,:));
%         y = X(j,:) - mean(X(j,:));
%         Px = sum(x.^2);
%         Py = sum(y.^2);
%         Pxy = sum(x.*y);
%         P(i,j) = Pxy./(sqrt(Px.*Py));       
%     end
% end



