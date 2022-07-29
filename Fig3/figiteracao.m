iter = 0:1:7;
sum8usersaux = load('Curvaint8.mat');
sum5usersaux = load('Curvaint5.mat');
sum2usersaux = load('Curvaint2.mat');

sum8userszfaux = load('sum8usersintzf.mat');
sum5userszfaux = load('sum5usersintzf.mat');
sum2userszfaux = load('sum2usersintzf.mat');

sum8users = sum8usersaux.sumexp./100;
sum5users = sum5usersaux.sumexp./100;
sum2users = sum2usersaux.sumexp./100;
sum8userszf = sum(sum8userszfaux.sumTotal)/length(sum8userszfaux.sumTotal)./100;
sum5userszf = sum(sum5userszfaux.sumTotal)/length(sum5userszfaux.sumTotal)./100;
sum2userszf = sum(sum2userszfaux.sumTotal)/length(sum2userszfaux.sumTotal)./100;

figure
hold on
a = plot(iter,sum8users,'Color',[0 0.4470 0.7410]);
b = plot(iter,sum5users,'Color',[0.4940 0.1840 0.5560]);
c = plot(iter,sum2users,'Color',[0.6350 0.0780 0.1840]);
d = plot(iter,sum8userszf,'Color',[0 0.4470 0.7410]);
e = plot(iter,sum5userszf,'-o','Color',[0.4940 0.1840 0.5560]);
f = plot(iter,sum2userszf,'-o','Color',[0.6350 0.0780 0.1840]);
a1 = plot(iter,sum8users,'s','Color',[0 0.4470 0.7410]);
b1 = plot(iter,sum5users,'s','Color',[0.4940 0.1840 0.5560]);
c1 = plot(iter,sum2users,'s','Color',[0.6350 0.0780 0.1840]);
d1 = plot(iter,sum8userszf,'o','Color',[0 0.4470 0.7410]);
a2 = plot(iter,-1*ones(1,8),'ks');
a3 = plot(iter,-1*ones(1,8),'ko');

h = [a2,a3,a,b,c];
 legend(h,'JPUA','ZFUA','$K$ = 8','$K$ = 5','$K$ = 2')
