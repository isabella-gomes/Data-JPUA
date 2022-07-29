%Definir quantidade de usuários em VLC e RF inicialmente
%Definir a alocação inicial
clear all
tic
elements = {0:1}; %cell array with N vectors to combine
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
numExp=50;
%IndTotal = zeros(numExp,7);
%CSTotal=zeros(numExp,7);
%sumexp = zeros(numExp,1);
%alocUsers = resultaux(aloc,:);
parfor q=1:numExp
[gammaSrf,noiseRF,nUsers,nLeds,Ar,dv,Psi_k,phi,Psic,phi12,l,Ts,r,gamma,eta,echarge,I_DC,B,Xamb,i_amp,gVLC,gEveVLC,gEveVLCtil,noiseEveVLC,nAntenas,hbRF,heRF,heRFtil]=CSSNR;
aloc = randi([1 length(resultaux)],1,1);
alocUsers = resultaux(aloc,:);
%alocUsers = randi([0,1],1,nUsers);
convergence=0;
iteracao =0;
UsersVLC=0;
UsersRF=0;
CSUser=zeros(1,nUsers);
    CLUser=zeros(1,nUsers);
    CEUser=zeros(1,nUsers);
    CSVLC=zeros(1,length(CSUser(alocUsers==1)));
    CSRF=zeros(1,length(CSUser(alocUsers==0)));
    hUsers=zeros(nLeds,nUsers);
    hVLC=zeros(nLeds,length(CSUser(alocUsers==1)));
    hRF=zeros(nAntenas,length(CSUser(alocUsers==0)));
    j=1;
    k=1;
for i=1:nUsers
    if UsersVLC<nLeds && alocUsers(i)==1
        hUsers(:,i)=gVLC(:,i);
        hVLC(:,j)=gVLC(:,i);
        UsersVLC=UsersVLC+1;
        j=j+1;
    elseif UsersRF<nAntenas
        hUsers(:,i)=hbRF(:,i);
        hRF(:,k)=hbRF(:,i);
        UsersRF=UsersRF+1;
        k=k+1;
    end
end
%noiseVLC = real(2.*gamma.*echarge.*eta.*hVLC'.*I_DC.*B + 4.*pi.*echarge.*Ar.*gamma.*Xamb.*(1-cos(Psic)).*B+i_amp^2.*B);
noiseVLC = 10^(-21);
noise = mean(noiseVLC./(gamma.*eta),'all');
noiseEve = mean(noiseEveVLC./(gamma.*eta),'all');
alphak = (2^(2*log(1-(-1)))./(2*pi*exp(1).*noise));
alphae = (2^(2*log(1-(-1)))./(2*pi*exp(1).*noiseEve));
betak = (1./(3*noise));
betae = (1./(3*noiseEve));

if UsersVLC ~= 0
[Wotimo]=ZF(nLeds,UsersVLC,hVLC,gEveVLC,I_DC);
end
if UsersRF ~= 0
[Votimo]=ZFRF(nAntenas,UsersRF,hRF,heRF,gammaSrf);
end
j=1;
k=1;

for i=1:nUsers
    if alocUsers(i)==1
     W = Wotimo;
     W(:,j)=[];
      CLUser(i) = abs(log2((sum(alphak.*((hVLC(:,j))'*Wotimo).^2)+1)/(sum(betak.*(((hVLC(:,j))'*W).^2))+ 1))/2);
      CEUser(i) = abs(log2((sum(alphae.*((gEveVLC)'*Wotimo).^2)+1)/(sum(betae.*(((gEveVLC)'*W).^2))+ 1))/2);
      CSVLC(j)=max(CLUser(i)-CEUser(i),0);
      j=1+j;
    elseif alocUsers(i)==0
        V = Votimo;
        V(:,k)=[];
        CLUser(i) = log2(1+(((hRF(:,k)'*Votimo(:,k)).^2))./(sum((((hRF(:,k)'*V)).^2))+ noiseRF));
        CEUser(i) = log2(1+(((heRF'*Votimo(:,k)).^2))./(sum(((heRF'*V).^2))+noiseRF));
        CSRF(k)=max(CLUser(i)-CEUser(i),0);
        k=1+k;
    end
      CSUser(i) = max(CLUser(i)-CEUser(i),0);
end
CSVLCi = CSVLC;
CSRFi = CSRF;
hUsersi = hUsers;
hVLCi = hVLC;
hRFi = hRF;
sumRF = sum(CSRF);
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
    if alocUsers(i)==0 && UsersVLC<nLeds
        CSaux = CSUseri;
        hRFaux = hRFi;
        hUsersi(:,i)=gVLC(:,i);
        UsersVLC=UsersVLC+1;
        hVLCi(:,UsersVLC)=gVLC(:,i);
       % noiseVLC = real(2.*gamma.*echarge.*eta.*hVLCi'.*I_DC.*B + 4.*pi.*echarge.*Ar.*gamma.*Xamb.*(1-cos(Psic)).*B+i_amp^2.*B);
        noisei = mean(noiseVLC./(gamma.*eta),'all');
        alphaki = (2^(2*log(1-(-1)))./(2*pi*exp(1).*noisei));
        betaki = (1./(3*noisei));
        [Wotimoi]=ZF(nLeds,UsersVLC,hVLCi,gEveVLC,I_DC);
        CLUseriVLC=zeros(1,UsersVLC);
        CEUseriVLC=zeros(1,UsersVLC);
        CSVLCaux = zeros(1,UsersVLC);
            for z=1:UsersVLC
                W = Wotimoi;
                W(:,z)=[];
                CLUseriVLC(z) = abs(log2((sum(alphaki.*((hVLCi(:,z))'*Wotimoi).^2)+1)/(sum(betaki.*(((hVLCi(:,z))'*W).^2))+ 1))/2);
                CEUseriVLC(z) = abs(log2((sum(alphae.*((gEveVLC)'*Wotimoi).^2)+1)/(sum(betae.*(((gEveVLC)'*W).^2))+ 1))/2);
                CSVLCaux(z)=max(CLUseriVLC(z)-CEUseriVLC(z),0);
            end
            for b = 1:length(hRFi)
            if hRFi(:,b)==hbRF(:,i)
                hRFi(:,b)=[];
                break
            end
           end
        %hRFi(:,j)=[];
        UsersRF=UsersRF-1;
        if UsersRF~=0
            [Votimoaux]=ZFRF(nAntenas,UsersRF,hRFi,heRF,gammaSrf);
            CLUseriRF=zeros(1,UsersRF);
            CEUseriRF=zeros(1,UsersRF);
            CSRFaux=zeros(1,UsersRF);
            for z=1:UsersRF
                Vaux = Votimoaux;
                Vaux(:,z)=[];
                CLUseriRF(z) = log2(1+(((hRFi(:,z)'*Votimoaux(:,z)).^2))./(sum((((hRFi(:,z)'*Vaux)).^2))+ noiseRF));
                CEUseriRF(z) = log2(1+(((heRF'*Votimoaux(:,z)).^2))./(sum(((heRF'*Vaux).^2))+noiseRF));
                CSRFaux(z)=max(CLUseriRF(z)-CEUseriRF(z),0);
            end
            sumRFi = sum(CSRFaux);
        else
            Votimoaux=0;
            CSRFaux = 0;
            sumRFi = 0;
        end
        if sumRFi + sum(CSVLCaux) >sumRF+ sumVLC       
           alocUsers(i)=1;
           %UsersRF=UsersRF-1;
           sumVLC=sum(CSVLCaux);
           sumRF=sumRFi;
           CSVLCi = CSVLCaux;
           CSRFi = CSRFaux;
           [CSUseri]=CalculoCS(alocUsers,hRFi,heRF,hVLCi,gEveVLC,Wotimoi,Votimoaux,noiseRF,betaki,betae,alphaki,alphae,nUsers);
        else 
            CSUseri = CSaux;
           hUsersi(:,i)=hUsers(:,i);
           hVLCi(:,UsersVLC)=[];
           hRFi = hRFaux;
           UsersVLC=UsersVLC-1;
           UsersRF=UsersRF+1;
           alphaki=alphak;
           betaki=betak;
           j=j+1;
        end
    elseif alocUsers(i)==1 && UsersRF<nAntenas
        CSaux = CSUseri;
        hVLCaux = hVLCi;
        hUsersi(:,i)=hbRF(:,i);
        UsersRF=UsersRF+1;
        hRFi(:,UsersRF)=hbRF(:,i);
        %[Vi] = ValorInicial(noiseRF,gammak,hRFi,heRF,UsersRF,nAntenas);
        [Votimoi]=ZFRF(nAntenas,UsersRF,hRFi,heRF,gammaSrf);
        CLUseriRF=zeros(1,UsersRF);
        CEUseriRF=zeros(1,UsersRF);
        CSRFaux=zeros(1,UsersRF);
            for z=1:UsersRF
                Vaux = Votimoi;
                Vaux(:,z)=[];
                CLUseriRF(z) = log2(1+(((hRFi(:,z)'*Votimoi(:,z)).^2))./(sum((((hRFi(:,z)'*Vaux)).^2))+ noiseRF));
                CEUseriRF(z) = log2(1+(((heRF'*Votimoi(:,z)).^2))./(sum(((heRF'*Vaux).^2))+noiseRF));
                CSRFaux(z)=max(CLUseriRF(z)-CEUseriRF(z),0);
            end
            for b = 1:length(hVLCi)
            if hVLCi(:,b)==gVLC(:,i)
                hVLCi(:,b)=[];
                break
            end
            end
        UsersVLC=UsersVLC-1;
        if UsersVLC~=0
            %noiseVLC = real(2.*gamma.*echarge.*eta.*hVLCi'.*I_DC.*B + 4.*pi.*echarge.*Ar.*gamma.*Xamb.*(1-cos(Psic)).*B+i_amp^2.*B);
            noisei = mean(noiseVLC./(gamma.*eta),'all');
            alphaki = (2^(2*log(1-(-1)))./(2*pi*exp(1).*noisei));
            betaki = (1./(3*noisei));
            [Wotimoaux]=ZF(nLeds,UsersVLC,hVLCi,gEveVLC,I_DC);
            CLUseriVLC = zeros(1,UsersVLC);
            CEUseriVLC = zeros(1,UsersVLC);
            CSVLCaux = zeros(1,UsersVLC);
            for f=1:UsersVLC
                Waux = Wotimoaux;
                Waux(:,f)=[];
                CLUseriVLC(f) = abs(log2((sum(alphaki.*((hVLCi(:,f))'*Wotimoaux).^2)+1)/(sum(betaki.*(((hVLCi(:,f))'*Waux).^2))+ 1))/21);
                CEUseriVLC(f) = abs(log2((sum(alphae.*((gEveVLC)'*Wotimoaux).^2)+1)/(sum(betae.*(((gEveVLC)'*Waux).^2))+ 1))/2);
                CSVLCaux(f)=max(CLUseriVLC(f)-CEUseriVLC(f),0);
            end
            sumVLCi = sum(CSVLCaux);
        else
             noisei = 0;
            alphaki = (2^(2*log(1-(-1)))./(2*pi*exp(1).*noisei));
            betaki = (1./(3*noisei));
            Wotimoaux=0;
            CSVLCaux = 0;
            sumVLCi = 0;
        end
        if  sumVLCi + sum(CSRFaux) > sumVLC + sumRF 
            alocUsers(i)=0;
            [CSUseri]=CalculoCS(alocUsers,hRFi,heRF,hVLCi,gEveVLC,Wotimoaux,Votimoi,noiseRF,betaki,betae,alphaki,alphae,nUsers);
                %UsersVLC=UsersVLC-1;
            sumRF=sum(CSRFaux);
            sumVLC=sumVLCi;
            CSRFi = CSRFaux;
            CSVLCi = CSVLCaux;
        else 
            hUsersi(:,i)=hUsers(:,i);
            hRFi(:,UsersRF)=[];
            hVLCi=hVLCaux;
            UsersRF=UsersRF-1;
            UsersVLC=UsersVLC+1;
            k=k+1;
            CSUseri=CSaux;
        end
    elseif alocUsers(i)==1 && UsersRF==nAntenas
        CSaux = CSUseri;
        hRFaux = hRFi;
        hUsersi(:,i)=hbRF(:,i);
        minRF = max(CSUseri(alocUsers==0));
        aux = 1;
        for a = 1:nUsers
            if alocUsers(a)==0 && CSUseri(a) == minRF
                ind =aux;
                break
            else 
                aux = aux+1;
            end
        end
        hRFi(hRFi==hbRF(:,ind))=[];
        hRFi=reshape(hRFi,nAntenas,UsersRF-1);
        hUsersi(:,ind) = gVLC(:,ind);
        for b = 1:length(hVLCi)
            if hVLCi(:,b)==gVLC(:,i)
                hVLCi(:,b)=gVLC(:,ind);
                break
            end
        end
        %hVLCi(hVLCi==gVLC(:,i))=gVLC(:,ind);
        %hVLCi = reshape(hVLCi,nLeds,UsersVLC);
        hRFi(:,UsersRF)=hbRF(:,i);
        %noiseVLC = real(2.*gamma.*echarge.*eta.*hVLCi'.*I_DC.*B + 4.*pi.*echarge.*Ar.*gamma.*Xamb.*(1-cos(Psic)).*B+i_amp^2.*B);
        noisei = mean(noiseVLC./(gamma.*eta),'all');
        alphaki = (2^(2*log(1-(-1)))./(2*pi*exp(1).*noisei));
        betaki = (1./(3*noisei));
        [Wotimoi]=ZF(nLeds,UsersVLC,hVLCi,gEveVLC,I_DC);
        CLUseriVLC=zeros(1,UsersVLC);
        CEUseriVLC=zeros(1,UsersVLC);
        CSVLCaux = zeros(1,UsersVLC);
            for z=1:UsersVLC
                W = Wotimoi;
                W(:,z)=[];
                CLUseriVLC(z) = abs(log2((sum(alphaki.*((hVLCi(:,z))'*Wotimoi).^2)+1)/(sum(betaki.*(((hVLCi(:,z))'*W).^2))+ 1))/2);
                CEUseriVLC(z) = abs(log2((sum(alphae.*((gEveVLC)'*Wotimoi).^2)+1)/(sum(betae.*(((gEveVLC)'*W).^2))+ 1))/2);
                CSVLCaux(z)=max(CLUseriVLC(z)-CEUseriVLC(z),0);
            end
        [Votimoaux]=ZFRF(nAntenas,UsersRF,hRFi,heRF,gammaSrf);
        CLUseriRF=zeros(1,UsersRF);
        CEUseriRF=zeros(1,UsersRF);
        CSRFaux=zeros(1,UsersRF);
            for z=1:UsersRF
                Vaux = Votimoaux;
                Vaux(:,z)=[];
                CLUseriRF(z) = log2(1+(((hRFi(:,z)'*Votimoaux(:,z)).^2))./(sum((((hRFi(:,z)'*Vaux)).^2))+ noiseRF));
                CEUseriRF(z) = log2(1+(((heRF'*Votimoaux(:,z)).^2))./(sum(((heRF'*Vaux).^2))+noiseRF));
                CSRFaux(z)=max(CLUseriRF(z)-CEUseriRF(z),0);
            end
        if sum(CSRFaux) + sum(CSVLCaux) >sumRF+ sumVLC       
           alocUsers(i)=0;
           alocUsers(ind)=1;
           [CSUseri]=CalculoCS(alocUsers,hRFi,heRF,hVLCi,gEveVLC,Wotimoi,Votimoaux,noiseRF,betaki,betae,alphaki,alphae,nUsers);
                %UsersRF=UsersRF-1;
           CSRFi = CSRFaux;
           CSVLCi = CSVLCaux;
           sumVLC=sum(CSVLCi);
           sumRF=sum(CSRFi);
        else 
           hUsersi(:,i)=gVLC(:,i);
           hUsersi(:,ind)=hbRF(:,ind);
           hRFi(:,UsersRF)=hbRF(:,ind);
           for b = 1:length(hVLCi)
            if hVLCi(:,b)==gVLC(:,ind)
                hVLCi(:,b)=gVLC(:,i);
                break
            end
           end
           alphaki=alphak;
           betaki=betak;
           CSUseri=CSaux;
        end
        elseif alocUsers(i)==0 && UsersVLC==nLeds
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
        hUsersi(:,ind) = hbRF(:,ind);
        for b = 1:length(hRFi)
            if hRFi(:,b)==hbRF(:,i)
                hRFi(:,b)=hbRF(:,ind);
                break
            end
        end
        %noiseVLC = real(2.*gamma.*echarge.*eta.*hVLCi'.*I_DC.*B + 4.*pi.*echarge.*Ar.*gamma.*Xamb.*(1-cos(Psic)).*B+i_amp^2.*B);
        noisei = mean(noiseVLC./(gamma.*eta),'all');
        alphaki = (2^(2*log(1-(-1)))./(2*pi*exp(1).*noisei));
        betaki = (1./(3*noisei));
        [Wotimoi]=ZF(nLeds,UsersVLC,hVLCi,gEveVLC,I_DC);
        CLUseriVLC=zeros(1,UsersVLC);
        CEUseriVLC=zeros(1,UsersVLC);
        CSVLCaux = zeros(1,UsersVLC);
            for z=1:UsersVLC
                W = Wotimoi;
                W(:,z)=[];
                CLUseriVLC(z) = abs(log2((sum(alphaki.*((hVLCi(:,z))'*Wotimoi).^2)+1)/(sum(betaki.*(((hVLCi(:,z))'*W).^2))+ 1))/2);
                CEUseriVLC(z) = abs(log2((sum(alphae.*((gEveVLC)'*Wotimoi).^2)+1)/(sum(betae.*(((gEveVLC)'*W).^2))+ 1))/2);
                CSVLCaux(z)=max(CLUseriVLC(z)-CEUseriVLC(z),0);
            end
        [Votimoaux]=ZFRF(nAntenas,UsersRF,hRFi,heRF,gammaSrf);
        CLUseriRF=zeros(1,UsersRF);
        CEUseriRF=zeros(1,UsersRF);
        CSRFaux=zeros(1,UsersRF);
            for z=1:UsersRF
                Vaux = Votimoaux;
                Vaux(:,z)=[];
                CLUseriRF(z) = log2(1+(((hRFi(:,z)'*Votimoaux(:,z)).^2))./(sum((((hRFi(:,z)'*Vaux)).^2))+ noiseRF));
                CEUseriRF(z) = log2(1+(((heRF'*Votimoaux(:,z)).^2))./(sum(((heRF'*Vaux).^2))+noiseRF));
                CSRFaux(z)=max(CLUseriRF(z)-CEUseriRF(z),0);
            end
        if sum(CSRFaux) + sum(CSVLCaux) >sumRF+ sumVLC       
           alocUsers(i)=1;
           alocUsers(ind)=0;
           [CSUseri]=CalculoCS(alocUsers,hRFi,heRF,hVLCi,gEveVLC,Wotimoi,Votimoaux,noiseRF,betaki,betae,alphaki,alphae,nUsers);
                %UsersRF=UsersRF-1;
           CSVLCi = CSVLCaux;
           CSRFi = CSRFaux;
           sumVLC=sum(CSVLCi);
           sumRF=sum(CSRFi);
        else 
           hUsersi(:,i)=hbRF(:,i);
           hUsersi(:,ind)=gVLC(:,ind);
           for b = 1:length(hVLCi)
            if hVLCi(:,b)==gVLC(:,i)
                hVLCi(:,b)=gVLC(:,ind);
                break
            end
           end
           for b = 1:length(hRFi)
            if hRFi(:,b)==hbRF(:,ind)
                hRFi(:,b)=hbRF(:,i);
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
    if UsersVLC>0
    [Wotimoi]=ZF(nLeds,UsersVLC,hVLCi,gEveVLC,I_DC);
    else
        Wotimo=0;
    end
    if UsersRF>0
    [Votimoaux]=ZFRF(nAntenas,UsersRF,hRFi,heRF,gammaSrf);
    else
        Votimoaux=0;
    end
    [CSUserfinal]=CalculoCS(alocUsers,hRFi,heRFtil,hVLCi,gEveVLCtil,Wotimoi,Votimoaux,noiseRF,betak,betae,alphak,alphae,nUsers);
    SumCS = sum(sumVLC)+sum(sumRF);
    sumcomp = sum(CSUserfinal);
else
    hUsers=hUsersi;
    CSRF=CSRFi;
    CSVLC=CSVLCi;
    CSUser=CSUseri;
end
end
%sumcompara(p) = sumcomp;
IndTotal(q,:)=alocUsers;
CSTotal(q,:) = CSUserfinal; 
sumexp(q,:) = sumcomp;
end
%save sumAlg2users(q).mat sumTotal
%save IndAlg2users(q).mat IndTotal

%save sumAlg7users4.mat sumcomp
%save CSTotalAlg7users4.mat CSUseri
%save IndTotalAlg7users3.mat alocUsers

save IndzfAlg1users1.mat IndTotal
save CSTotalzfAlg1users1.mat CSTotal
save sumzfAlg1users1.mat sumexp
toc
