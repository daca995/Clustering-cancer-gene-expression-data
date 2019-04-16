function [IDX]= fmt_NG_spec_clust(dist,k,sigma)
% k je broj klastera 
%Y=[]; 
[m,n]=size(dist); 


% nekad se sigma uzima kao mean od matrice dist, mozemo razmotiri i tu
% kombinaciju
%calculo do sigma, depois de fazer clustering �s linhas de Y, escolher o valor de modo a minimizar a distor��o 
%(tightest (smallest distortion) clusters) -> a fazer


A=exp(-dist.^2/(2.*sigma^2));

A(sub2ind([m n],[1:m],[1:m]))=0;

th=1e-10;
lixo=length(find(sum(A)<th));
if lixo/n>.10 
   fprintf(1,'Mais de %2.1f amostras (%i) sao zero ou proximo disso!\n',lixo/n,lixo);
	fprintf(1,'A sair...\n');
    IDX=[]
   return   
end


D=diag(sum(A,2));

aux=D^(-0.5);

L=aux*A*aux;

[vect,diag_eig] = eig(L);
%[V,D,FLAG] = eigs(A,B,K,SIGMA,OPTIONS)

%estou a considerar o uso dos squared eigenvalues... a fazer
[sorted,indices] = sort(diag_eig(sub2ind([m n],[1:m],[1:m])));
%ordem ascendente (end � o maior)

i=0;
j=0;
eig_usados=[];
X=[];


while i<=k & j<=m
   X(1:m,i+1)=vect(:,indices(end-j));
   eig_usados(i+1)=diag_eig(indices(end-j),indices(end-j));
   %verificar se s�o iguais:
   if i>0
	   if sum(sum(ones([1 i]).*diag_eig(indices(end-j),indices(end-j))==eig_usados(1:i) ,2)) == 0
         i=i+1;
       end
    else
       i=i+1;
   end
    
   j=j+1;
end


%so extendi X para observar o eigen gap
X=X(1:m,1:k);
eig_usados=eig_usados(1:end-1);




%normalizar as linhas
Y=diag(1./sqrt(sum(X.^2,2)))*X;

%comentei isto para acelerar o processo
%plotmatrix(Y);

%print(gcf,'-depsc','-noui',[fich(1:end-4) '-A-' num2str(k) '-' num2str(sigma) '.eps']);
%cfig=gcf;
%delete(cfig);

%tratar cada linha de Y como sendo uma amostra no espa�o R^k e agrupa-las em k clusters usando kmedias

%k_medias_with_seed_vns.m  --- faz clustering em dois passos
%k_medias_with_seed_vnsc --- faz clustering em dois passos e retorna tb os centroides


%we let the first centroid be randomly chosen row of Y, and then repeatedly choose as the
%next centroid the row of Y that is closest to being 90� from all the centroids already picked.
%seed_order=(1:k);

seed_order(1)=ceil(rand*m);

for i=2:k
   prod=abs(Y*Y(seed_order(1:i-1),:)');
   %REPMAT Replicate and tile an array
   %prod=abs(repmat(Y,[1 i-1])*Y(seed_order(1:i-1),:)');
   
   %prod � uma matriz que vai representar o produto interno de todos os elementos de y com os centroides ja escolhidos
   %o proximo a ser escolhido � o estiver mais perto de estar a 90� dos centroides anteriores, ou seja,
   %o que tiver valor minimo de (soma de) produtos internos (1� coluna produto interno,pi, de Y com o 1� centroide,
   %2� pi de Y c o 2� centroide
   
   %[sorted,indices]=min(sum(prod,2)); 
   %mudei para o sort para ter a possibilidade de escolher outro que nao o minimo
   %20 Maio
   [sorted,indices]=sort(sum(prod,2));
   
   %aux=find(indices~=seed_order(1:i-1));
   j=1;
   while j<=m
      ok=1;
      %ver se j� foi escolhido anteriormente
      for h=1:i-1
         if indices(j)==seed_order(h)
				ok=0;
         end
      end
      if ok==0 
         %se sim entao escolhemos o pr�ximo mais perpendicular
         j=j+1;
      else
         seed_order(i)=indices(j);
         break
      end
   end
end

%they used kmeans with only 1 run (no restart) 
%[nsamples_in_cluster,clusters_m]=k_medias_with_seed_vns(Y,k,seed_order);

%[ns_in_cl,clusters_m, centroides]=k_medias_with_seed_vnsc(Y,k,seed_order);
warning('OFF','stats:kmeans:EmptyCluster');
warning('OFF','stats:kmeans:FailedToConverge');
opts = statset('Display','off','MaxIter',1);%,'Display','final')
IDX = kmeans(Y,k,'start',Y(seed_order,:),'distance','sqEuclidean','EmptyAction','drop','Options',opts);

%% iscrtavanje

% gscatter(Y(:,1),Y(:,2),IDX,'br','xo'), title('gscatter')
% figure, gplot(A,Y), title('gplot')
% PlotGraph2(A, IDX), title('PlotGraph2')

function D = eucldistm(A,B)
%EUCLDISTM  Euclidean distance matrix
%
%        D = EUCLDISTM(A,B)
%
% A specialized function for computing the squared Euclidean
% distance D between datasets A and B. This is mainly for
% computational speed, so it is light-weight, without any checking.
% Normal users will probably use myproxm.
%
% Mainly for internal use.
%
% see also: myproxm

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org; Ana Fred
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

D = sqrt(sqeucldistm(A,B));

return

function D = sqeucldistm(A,B)
%SQEUCLDISTM Square Euclidean distance matrix
%
%        D = SQEUCLDISTM(A,B)
%
% A specialized function for computing the squared Euclidean
% distance D between datasets A and B. This is mainly for
% computational speed, so it is light-weight, without any checking.
% Normal users will probably use myproxm.
%
% Mainly for internal use.
%
% see also: myproxm

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

nra = size(A,1);
nrb = size(B,1);
dA = sum(A.*A,2);
dB = sum(B.*B,2);
D = repmat(dA,1,nrb) + repmat(dB',nra,1) -2*A*B';
% no negative distances:
D(D<0) = 0;

return
