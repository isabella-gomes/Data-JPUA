%Definir quantidade de usuários em VLC e RF inicialmente
%Definir a alocação inicial
clear all
tic
elements = {0:1, 0:1, 0:1, 0:1, 0:1,0:1, 0:1, 0:1}; %cell array with N vectors to combine
 combinations = cell(1, numel(elements)); %set up the varargout result
 [combinations{:}] = ndgrid(elements{:});
 combinations = cellfun(@(x) x(:), combinations,'uniformoutput',false); %there may be a better way to do this
 result = [combinations{:}]; % NumberOfCombinations by N matrix. Each row is unique.
 resultaux = result;
 for i=length(result):-1:1
     if length(result(result(i,:)==1))>4
         resultaux(i,:)=[];
     elseif length(result(result(i,:)==0))>4
         resultaux(i,:)=[];
     end
 end
rate = 1;
numExp=30;
%IndTotal = zeros(numExp,7);
%CSTotal=zeros(numExp,7);
%sumexp = zeros(numExp,1);
%aloc = randi([1 length(resultaux)],1,1);
%alocUsers = resultaux(aloc,:);
parfor q=1:numExp
[gammaSrf,noiseRF,nUsers,nLeds,Ar,dv,Psi_k,phi,Psic,phi12,l,Ts,r,gamma,eta,echarge,I_DC,B,Xamb,i_amp,gVLC,gEveVLC,gEveVLCtil,noiseEveVLC,nAntenas,hbRF,heRF,heRFtil]=CSSNR;
%aloc = perms([1 0 0 0 0]);
%indice = randi(length(aloc));
%alocUsers = aloc(indice,:);
alocUsers=1;
%alocUsers = randi([0,1],1,nUsers);
convergence=0;
iteracao =0;
UsersVLC=0;
UsersRF=0;
CSUser=zeros(1,nUsers);
    CLUser=zeros(1,nUsers);
    CEUser=zeros(1,nUsers);
    CSVLC=zeros(1,length(CSUser(alocUsers==1)));
    %CSVLC=0;
    %CSRF=0;
    hUsers=zeros(nLeds,nUsers);
    hVLC=zeros(nLeds,length(CSUser(alocUsers==1)));
    %hVLC = zeros(nLeds,1);
    %hRF = zeros(nAntenas,1);
    j=1;
for i=1:nUsers
    if UsersVLC<nLeds && alocUsers(i)==1
        hUsers(:,i)=gVLC(:,i);
        hVLC(:,j)=gVLC(:,i);
        UsersVLC=UsersVLC+1;
        j=j+1;
    end
end
noiseVLC = 10^(-21);
alphak = (2^(2*log(1-(-1)))./(2*pi*exp(1).*noiseVLC));
alphae = (2^(2*log(1-(-1)))./(2*pi*exp(1).*noiseEveVLC));
betak = (1./(3*noiseVLC));
betae = (1./(3*noiseEveVLC));

if UsersVLC ~= 0
[Wotimo]=algoritmoVLC(noiseVLC,noiseEveVLC,rate,I_DC,gamma,eta,hVLC,gEveVLC,UsersVLC,nLeds);
end
j=1;
k=1;

for i=1:nUsers
    if alocUsers(i)==1
     W = Wotimo;
     W(:,j)=[];
      CLUser(i) = log2((sum(alphak.*((hVLC(:,j))'*Wotimo).^2)+1)/(sum(betak.*(((hVLC(:,j))'*W).^2))+ 1))/2;
      CEUser(i) = log2((sum(alphae.*((gEveVLC)'*Wotimo).^2)+1)/(sum(betae.*(((gEveVLC)'*W).^2))+ 1))/2;
      CSVLC(j)=max(CLUser(i)-CEUser(i),0);
      j=1+j;
    elseif alocUsers(i)==0
        CLUser(i) = 0;
        CEUser(i) = 0;
    end
      CSUser(i) = max(CLUser(i)-CEUser(i),0);
end
CSVLCi = CSVLC;
hUsersi = hUsers;
hVLCi = hVLC;
sumVLC = sum(CSVLC);
%% Iteração 
while convergence==0
    CSUseri = CSUser;
    alocTest=alocUsers;
    CLUseri=zeros(1,nUsers);
    CEUseri=zeros(1,nUsers);
    j=1;
    k=1;
%% Calculo do payoff
for i = 1:nUsers
        if alocUsers(i)==0 && UsersVLC==nLeds
        CSaux = CSUseri;
        hVLCaux = hVLCi;
        hUsersi(:,i)=gVLC(:,i);
        minVLC = max(CSUseri(alocUsers==1));
        aux = 1;
        for a = 1:nUsers
            if alocUsers(a)==1 && CSUseri(a) == minVLC
                ind =aux;
                break
            else 
                aux = aux+1;
            end
        end
        %ind = aux;
        %ind = find(CSUseri==minVLC,1);
        for b = 1:length(hVLCi)
            if hVLCi(:,b)==gVLC(:,ind)
                hVLCi(:,b)=gVLC(:,i);
                break
            end
        end
        noiseVLCi = 10^(-21);
        alphaki = (2^(2*log(1-(-1)))./(2*pi*exp(1).*noiseVLC));
        betaki = (1./(3*noiseVLC));
        [Wotimoi]=algoritmoVLC(noiseVLCi,noiseEveVLC,rate,I_DC,gamma,eta,hVLCi,gEveVLC,UsersVLC,nLeds);
        CLUseriVLC=zeros(1,UsersVLC);
        CEUseriVLC=zeros(1,UsersVLC);
        CSVLCaux = zeros(1,UsersVLC);
            for z=1:UsersVLC
                W = Wotimoi;
                W(:,z)=[];
                CLUseriVLC(z) = log2((sum(alphaki.*((hVLCi(:,z))'*Wotimoi).^2)+1)/(sum(betaki.*(((hVLCi(:,z))'*W).^2))+ 1))/2;
                CEUseriVLC(z) = log2((sum(alphae.*((gEveVLC)'*Wotimoi).^2)+1)/(sum(betae.*(((gEveVLC)'*W).^2))+ 1))/2;
                CSVLCaux(z)=max(CLUseriVLC(z)-CEUseriVLC(z),0);
            end
        if sum(CSVLCaux)>sumVLC       
           alocUsers(i)=1;
           alocUsers(ind)=0;
           %[CSUseri]=CalculoCS(alocUsers,hRFi,heRF,hVLCi,gEveVLC,Wotimoi,Votimoaux,noiseRF,betaki,betae,alphaki,alphae,nUsers);
                %UsersRF=UsersRF-1;
           CSVLCi = CSVLCaux;
           sumVLC=sum(CSVLCi);
        else 
           for b = 1:length(hVLCi)
            if hVLCi(:,b)==gVLC(:,i)
                hVLCi(:,b)=gVLC(:,ind);
                break
            end
           end
           alphaki=alphak;
           betaki=betak;
           CSUseri=CSaux;
        end
        end
end
if alocUsers==alocTest
    convergence=1;
    [Wotimoi]=algoritmoVLC(noiseVLC,noiseEveVLC,rate,I_DC,gamma,eta,hVLCi,gEveVLC,UsersVLC,nLeds);
    CSVLCfim = zeros(1,UsersVLC);
    for z=1:UsersVLC
                W = Wotimoi;
                W(:,z)=[];
                CLUseriVLC(z) = log2((sum(alphak.*((hVLCi(:,z))'*Wotimoi).^2)+1)/(sum(betak.*(((hVLCi(:,z))'*W).^2))+ 1))/2;
                CEUseriVLC(z) = log2((sum(alphae.*((gEveVLCtil)'*Wotimoi).^2)+1)/(sum(betae.*(((gEveVLCtil)'*W).^2))+ 1))/2;
                CSVLCfim(z)=max(CLUseriVLC(z)-CEUseriVLC(z),0);
    end
    sumVLCfim = sum(CSVLCfim);
    SumCS = sum(sumVLCfim);
else
    %SumCS = 0;
    hUsers=hUsersi;
    CSVLC=CSVLCi;
end
end
%sumcompara(p) = sumcomp;
IndTotal(q,:)=alocUsers;
CSTotal(q,:) = CSVLCfim; 
sumexp(q,:) = sum(CSVLCfim);
end
%save sumAlg2users(q).mat sumTotal
%save IndAlg2users(q).mat IndTotal

%save sumAlg7users4.mat sumcomp
%save CSTotalAlg7users4.mat CSUseri
%save IndTotalAlg7users3.mat alocUsers

save IndAlg1usersVLC1.mat IndTotal
save CSTotalAlg1VLC1.mat CSTotal
save sumAlg1usersVLC1.mat sumexp
toc
