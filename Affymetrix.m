%PROJEKAT IZ PREPOZNAVANJA OBLIKA
% Klasterizacija gena
% Daniela Kotur EE95/2014


tic
name = {'armstrong_v1','armstrong_v2','bhattacharjee_2001','chowdary','dyrskjot','golub_v1','golub_v2','gordon','laiho','nutt_v1', 'nutt_v2','nutt_v3', 'pomeroy_v1', 'pomeroy_v2','ramaswamy','shipp','singh','su', 'west','yeoh_v1','yeoh_v2'};
protokol = {'Helingerovo_rastojanje','E_RangeR', 'E_RangeC', 'E_ZscoreR', 'E_ZscoreC', 'Pirsonova_korelacija'};

%Spektralno klasterovanje 

for z = 1:6 %6 protokola
   for j = 1:21 %21 set podataka
       ime = strcat(name{j},'.txt'); %niz koji sadrzi imena setova podataka
       p_ime = protokol{z};  %niz koji sadrzi imena protokola 
       n=name{j}; 
       f = str2func(p_ime); 
       M = importdata(ime);  
       dist = f(M); %uzima svaki poseban protokol(rastojanja, uz prethodno normalizovane ili nenormalizovane podatke)
       sigma = mean(dist(:)); %sigma se uzima kao srednja vrijednost rastojanja
       labele = importdata(strcat(name{j}, '_labels.txt'));
       k = max(labele); %broj klastera, odnosno klasa u svakom setu podataka
    
   % ari=[];
     for i = 1:15 %15 iteracija
         IDX = fmt_NG_spec_clust(dist, k, sigma); %spektralno klasterovanje
         ari(i) = rand_index(IDX, labele, 'adjusted'); %Adjusted rand index, mjera validacije
     end
       save(['ari_' n '_'  p_ime '.mat'],'ari'); %cuvanje rezultata spektralnog klasterovanja svih 15 iteracija svakog seta podataka
       
   end
end

%Iscrtavanje rezultata pomocu stubica

ARI=[];
for k=0:5:15 %radi predstavljanja rezultata kao 4 grupe po 5 setova podataka
   for j=k+1:k+5 %21 set podataka 
       for i= 1:6 %6 protokola
           p_ime = protokol{i};
           n = name{j};
           b = load(['ari_' n '_'  p_ime '.mat'],'ari'); %ucitavanje svih sacuvanih rezultata
           box = b.ari;  
           ARI = [ARI; box];
    
       end
   end
   
   ARI = ARI';
   m = mean(ARI); 
   w = [m(1) m(2) m(3) m(4) m(5) m(6); m(7) m(8) m(9) m(10) m(11) m(12); m(13) m(14) m(15) m(16) m(17) m(18); m(19) m(20) m(21) m(22) m(23) m(24); m(25) m(26) m(27) m(28) m(29) m(30)];
   bar(w); %iscrtavanje 
   ARI=[];
   
end

 
for i= 1:6 
    p_ime = protokol{i};
    n = name{21}; %posljednji set
    b= load(['ari_' n '_'  p_ime '.mat'],'ari');
    box = b.ari;
    ARI = [ARI; box];   
end

 ARI = ARI';
 m = mean(ARI);
 bar(m);
 ARI = [];

%Iscrtavanje rezultata pomocu boxplota

for i = 1:6
    for j= 1:21
        p_ime = protokol{i};
        n = name{j};
        b = load(['ari_' n '_'  p_ime '.mat'],'ari');
        box = b.ari;
        ARI = [ARI; box];
    end
        
    ARI = ARI';
    figure, boxplot(ARI, 'Labelorientation', 'inline','Labels', {'Armstrong_v1', 'Armstrong_v2', 'Bhattacharjee', 'Chowdary','Dyrskjot','Golub_v1', 'Golub_v2', 'Gordon', 'Laiho', 'Nutt_v1', 'Nutt_v2', 'Nutt_v3', 'Pomeroy_v1', 'Pomeroy_v2', 'Ramaswamy', 'Shipp', 'Singh', 'Su', 'West', 'Yeoh_v1', 'Yeoh_v2'} );
    title(protokol{i});
    ylabel('ARI');
    ARI = [];
     
end

toc

