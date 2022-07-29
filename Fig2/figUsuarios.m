numUsers = 1:1:8;
users8aux = load('sumAlg8users.mat');
%users7aux1 = load('sumAlg7users1.mat');
users7aux = load('sumAlg7users1.mat');
%users6aux3 = load('sumAlg6users3.mat');
%users6aux2 = load('sumAlg6users2.mat');
%users6aux1 = load('sumAlg6users1.mat');
users6aux = load('sumAlg6users1.mat');
%users5aux2 = load('sumAlg5users2.mat');
%users5aux1 = load('sumAlg5users1.mat');
users5aux = load('sumAlg5users1.mat');
%users4aux2 = load('sumAlg4users2.mat');
users4aux = load('sumAlg4users.mat');
%users4aux = load('sumAlg4users1.mat');
%users3aux2 = load('sumAlg3users2.mat');
%users3aux1 = load('sumAlg3users1.mat');
users3aux = load('sumAlg3users1.mat');
users2aux = load('sumAlg2users1.mat');
%users2aux = load('sumAlg2users.mat');
users1aux = load('sumAlg1users1.mat');
%users1aux = load('sumAlg1users.mat');


users8zfaux = load('sumzfAlg8users.mat');
users7zfaux = load('sumzfAlg7users1.mat');
users6zfaux = load('sumzfAlg6users1.mat');
users5zfaux = load('sumzfAlg5users1.mat');
users4zfaux = load('sumzfAlg4users1.mat');
users3zfaux = load('sumzfAlg3users1.mat');
users2zfaux = load('sumzfAlg2users1.mat');
users1zfaux = load('sumzfAlg1users1.mat');

users8randaux = load('sumAlg8usersrand.mat');
users7randaux = load('sumAlg7usersrand.mat');
users6randaux = load('sumAlg6usersrand.mat');
users5randaux = load('sumAlg5usersrand.mat');
users4randaux = load('sumAlg4usersrand.mat');
users3randaux = load('sumAlg3usersrand.mat');
users2randaux = load('sumAlg2usersrand.mat');
users1randaux = load('sumAlg1usersrand.mat');

users8vlcaux = load('sumAlg8usersVLC.mat');
users7vlcaux = load('sumAlg7usersVLC.mat');
users6vlcaux = load('sumAlg6usersVLC.mat');
users5vlcaux = load('sumAlg5usersVLC.mat');
users4vlcaux = load('sumAlg4usersVLC.mat');
users3vlcaux = load('sumAlg3usersVLC.mat');
users2vlcaux = load('sumAlg2usersVLC.mat');
users1vlcaux = load('sumAlg1usersVLC.mat');

users8rfaux = load('sumAlg8usersRF.mat');
users7rfaux = load('sumAlg7usersRF.mat');
users6rfaux = load('sumAlg6usersRF.mat');
users5rfaux = load('sumAlg5usersRF1.mat');
users4rfaux = load('sumAlg4usersRF1.mat');
users3rfaux = load('sumAlg3usersRF.mat');
users2rfaux = load('sumAlg2usersRF.mat');
users1rfaux = load('sumAlg1usersRF.mat');


users8vlc = (sum(users8vlcaux.sumexp))/(length(users8vlcaux.sumexp));
users7vlc = (sum(users7vlcaux.sumexp))/(length(users7vlcaux.sumexp));
users6vlc = (sum(users6vlcaux.sumexp))/(length(users6vlcaux.sumexp));
users5vlc = (sum(users5vlcaux.sumexp))/(length(users5vlcaux.sumexp));
users4vlc = (sum(users4vlcaux.sumexp))/(length(users4vlcaux.sumexp));
users3vlc = sum(users3vlcaux.sumexp)/(length(users3vlcaux.sumexp));
users2vlc = sum(users2vlcaux.sumexp)/length(users2vlcaux.sumexp);
users1vlc = sum(users1vlcaux.sumexp)/length(users1vlcaux.sumexp);

users8rf = (sum(users8rfaux.sumexp))/(length(users8rfaux.sumexp));
users7rf = (sum(users7rfaux.sumexp))/(length(users7rfaux.sumexp));
users6rf = (sum(users6rfaux.sumexp))/(length(users6rfaux.sumexp));
users5rf = (sum(users5rfaux.sumexp))/(length(users5rfaux.sumexp));
users4rf = (sum(users4rfaux.sumexp))/(length(users4rfaux.sumexp));
users3rf = sum(users3rfaux.sumexp)/(length(users3rfaux.sumexp));
users2rf = sum(users2rfaux.sumexp)/length(users2rfaux.sumexp);
users1rf = sum(users1rfaux.sumexp)/length(users1rfaux.sumexp);

users8 = (sum(users8aux.sumexp))/(length(users8aux.sumexp));
users7 = (sum(users7aux.sumexp)+sum(users7aux1.sumexp))/(length(users7aux.sumexp)+length(users7aux1.sumexp));
users6 = (sum(users6aux.sumexp))/(length(users6aux.sumexp));
users5 = (sum(users5aux.sumexp))/(length(users5aux.sumexp));
users4 = (sum(users4aux.sumexp))/(length(users4aux.sumexp));
users3 = (sum(users3aux.sumexp))/(length(users3aux.sumexp));
users2 = (sum(users2aux.sumexp))/(length(users2aux.sumexp));
users1 = (sum(users1aux.sumexp))/(length(users1aux.sumexp));

users8zf = sum(users8zfaux.sumexp)/length(users8zfaux.sumexp);
users7zf = sum(users7zfaux.sumexp)/length(users7zfaux.sumexp);
users6zf = sum(users6zfaux.sumexp)/length(users6zfaux.sumexp);
users5zf = sum(users5zfaux.sumexp)/length(users5zfaux.sumexp);
users4zf = sum(users4zfaux.sumexp)/length(users4zfaux.sumexp);
users3zf = sum(users3zfaux.sumexp)/length(users3zfaux.sumexp);
users2zf = sum(users2zfaux.sumexp)/length(users2zfaux.sumexp);
users1zf = sum(users1zfaux.sumexp)/length(users1zfaux.sumexp);

users8rand = sum(users8randaux.sumexp)/length(users8randaux.sumexp);
users7rand = sum(users7randaux.sumexp)/length(users7randaux.sumexp);
users6rand = sum(users6randaux.sumexp)/length(users6randaux.sumexp);
users5rand = (sum(users5randaux.sumexp))/(length(users5randaux.sumexp));
users4rand = sum(users4randaux.sumexp)/length(users4randaux.sumexp);
users3rand = sum(users3randaux.sumexp)/length(users3randaux.sumexp);
users2rand = sum(users2randaux.sumexp)/length(users2randaux.sumexp);
users1rand = sum(users1randaux.sumexp)/length(users1randaux.sumexp);

 usersrf = [users1rf users2rf users3rf users4rf users5rf users6rf users7rf users8rf]./100;
 userszf = [users1zf users2zf users3zf users4zf users5zf users6zf users7zf users8zf]./100;
 usersvlc = [users1vlc users2vlc users3vlc users4vlc users5vlc users6vlc users7vlc users8vlc]./100;
 usersrand = [users1rand users2rand users3rand users4rand users5rand users6rand users7rand users8rand]./100;

figure
users = [users1 users2 users3 users4 users5 users6 users7 users8]./100;
hold on
plot(numUsers,users,'Color',[0 0.4470 0.7410])
plot(numUsers,userszf,'Color',[0.4940 0.1840 0.5560])
plot(numUsers,usersrand,'Color',[0.4660 0.6740 0.1880])
plot(numUsers,usersvlc,'Color',[0.8500 0.3250 0.0980])
plot(numUsers,usersrf,'Color',[0.9290 0.6940 0.1250])
legend('JPUA', 'ZFUA', 'Random Association','VLC Stand-Alone', 'RF Stand-Alone')
%hold on
%semilogy(numUsers,userszf,'b-o')
%hold on
%semilogy(numUsers,usersrand,'r-o')
%hold on
%plot(numUsers,users,'k-o')
%hold on
%plot(numUsers,userszf,'b-o')
%plot(numUsers,usersrand,'r-o')
