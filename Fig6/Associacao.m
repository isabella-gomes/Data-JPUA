%Definir quantidade de usuários em VLC e RF inicialmente
%Definir a alocação inicial
clear all
tic
elements = {0:1,0:1,0:1,0:1,0:1}; %cell array with N vectors to combine
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
numExp=40;
%IndTotal = zeros(numExp,7);
%CSTotal=zeros(numExp,7);
%sumexp = zeros(numExp,1);
%alocUsers = resultaux(aloc,:);

for q=1:numExp
phi12 = pi/18:pi/18:(2*pi)/3; % angulo de meia intensidade do led
Sa = length(phi12);
%alocUsers = randi([0,1],1,nUsers);
parfor S=1:length(phi12)
users = -2.5 + (2.5+2.5)*rand(5,2); %Posicionando os usuários aleatóriamente
Eve = -2.5 + (2.5+2.5)*rand(1,2);
[Leds,gammaSrf,noiseRF,nUsers,nLeds,Ar,dv,Psi_k,phi,Psic,l,Ts,r,gamma,eta,echarge,I_DC,B,Xamb,i_amp,hbRF,heRF,heRFtil,nAntenas,Deltarf,Deltavlc]=CSSNR(phi12(S),users,Eve);
distsv = zeros(nLeds,nUsers);
distEvevlc = zeros(nLeds,1);
Raio = dv.*tan(phi12(S));
th = 0:pi/50:2*pi;
%Posicionando Cobertura de cada fonte
xLed1 = Raio.*cos(th)+Leds(1,1);
yLed1 = Raio.*sin(th)+Leds(1,2);
xLed2 = Raio.*cos(th)+Leds(2,1);
yLed2 = Raio.*sin(th)+Leds(2,2);
xLed3 = Raio.*cos(th)+Leds(3,1);
yLed3 = Raio.*sin(th)+Leds(3,2);
xLed4 = Raio.*cos(th)+Leds(4,1);
yLed4 = Raio.*sin(th)+Leds(4,2);
%% Canal VLC
%inpolygon - Definir se usuário esta interior ou na borda de cobertura do LED 
    [in1,on1] = inpolygon(users(:,1),users(:,2),xLed1,yLed1); %usuários legitimos
    [in2,on2] = inpolygon(users(:,1),users(:,2),xLed2,yLed2); %verificação circunferencia de cobertura
    [in3,on3] = inpolygon(users(:,1),users(:,2),xLed3,yLed3);
    [in4,on4] = inpolygon(users(:,1),users(:,2),xLed4,yLed4);
    for i=1:nUsers
    usersVLC = repmat(users(i,:),nLeds,1);
    dists =  sqrt(((usersVLC(:,1)-Leds(:,1)).^2)+((usersVLC(:,2)-Leds(:,2)).^2));
    distsv(:,i) = sqrt(dv^2 + dists.^2);
    end
    State1 = [in1 in2 in3 in4]';
    [ine1,one1] = inpolygon(Eve(1),Eve(2),xLed1,yLed1); %canal do espião
    [ine2,one2] = inpolygon(Eve(1),Eve(2),xLed2,yLed2);
    [ine3,one3] = inpolygon(Eve(1),Eve(2),xLed3,yLed3);
    [ine4,one4] = inpolygon(Eve(1),Eve(2),xLed4,yLed4);
    eveVLC = repmat(Eve,nLeds,1);
    distsEve =  sqrt(((eveVLC(:,1)-Leds(:,1)).^2)+((eveVLC(:,2)-Leds(:,2)).^2));
    distEvevlc(:,1) = sqrt(dv^2 + distsEve.^2);
    StateEve = [ine1;ine2;ine3;ine4];
gVLC= real(State1.*(Ar./(distsv.^2)).*((l+1)/(2*pi)).*cos(phi)^l.*Ts.*(r^2/(sin(Psic)^2)).*cos(Psi_k));
gEveVLC =real(StateEve.*(Ar./(distEvevlc.^2)).*((l+1)/(2*pi)).*cos(phi)^l.*Ts.*(r^2/(sin(Psic)^2)).*cos(Psi_k));
gEveVLCtil = gEveVLC + sqrt(Deltavlc).*randn(length(gEveVLC),1);
noiseEveVLC = 10^(-21);
aloc = randi([1 length(resultaux)],1,1);
alocUsers = resultaux(aloc,:);
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

noiseVLC = 10^(-21);
alphak = (2^(2*log(1-(-1)))./(2*pi*exp(1).*noiseVLC));
alphae = (2^(2*log(1-(-1)))./(2*pi*exp(1).*noiseEveVLC));
betak = (1./(3*noiseVLC));
betae = (1./(3*noiseEveVLC));

if UsersVLC ~= 0
[Wotimo]=algoritmoVLC(noiseVLC,noiseEveVLC,rate,I_DC,gamma,eta,hVLC,gEveVLC,UsersVLC,nLeds);
end
if UsersRF ~= 0
[Votimo]=algoritmoRFIterativo(noiseRF,rate,gammaSrf,hRF,heRF,UsersRF,nAntenas);
end
j=1;
k=1;

for i=1:nUsers
    if alocUsers(i)==1
     W = Wotimo;
     W(:,j)=[];
      CLUser(i) = (log2((sum(alphak.*((hVLC(:,j))'*Wotimo).^2)+1)/(sum(betak.*(((hVLC(:,j))'*W).^2))+ 1))/2);
      CEUser(i) = (log2((sum(alphae.*((gEveVLC)'*Wotimo).^2)+1)/(sum(betae.*(((gEveVLC)'*W).^2))+ 1))/2);
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
        noiseVLC = 10^(-21);
        alphaki = (2^(2*log(1-(-1)))./(2*pi*exp(1).*noiseVLC));
        betaki = (1./(3*noiseVLC));
        [Wotimoi]=algoritmoVLC(noiseVLC,noiseEveVLC,rate,I_DC,gamma,eta,hVLCi,gEveVLC,UsersVLC,nLeds);
        CLUseriVLC=zeros(1,UsersVLC);
        CEUseriVLC=zeros(1,UsersVLC);
        CSVLCaux = zeros(1,UsersVLC);
            for z=1:UsersVLC
                W = Wotimoi;
                W(:,z)=[];
                CLUseriVLC(z) = (log2((sum(alphaki.*((hVLCi(:,z))'*Wotimoi).^2)+1)/(sum(betaki.*(((hVLCi(:,z))'*W).^2))+ 1))/2);
                CEUseriVLC(z) = (log2((sum(alphae.*((gEveVLC)'*Wotimoi).^2)+1)/(sum(betae.*(((gEveVLC)'*W).^2))+ 1))/2);
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
            [Votimoaux]=algoritmoRFIterativo(noiseRF,rate,gammaSrf,hRFi,heRF,UsersRF,nAntenas);
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
        [Votimoaux]=algoritmoRFIterativo(noiseRF,rate,gammaSrf,hRFi,heRF,UsersRF,nAntenas);
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
            for b = 1:length(hVLCi)
            if hVLCi(:,b)==gVLC(:,i)
                hVLCi(:,b)=[];
                break
            end
            end
        UsersVLC=UsersVLC-1;
        if UsersVLC~=0
            alphaki = (2^(2*log(1-(-1)))./(2*pi*exp(1).*noiseVLC));
            betaki = (1./(3*noiseVLC));
            [Wotimoi]=algoritmoVLC(noiseVLC,noiseEveVLC,rate,I_DC,gamma,eta,hVLCi,gEveVLC,UsersVLC,nLeds);
            CLUseriVLC = zeros(1,UsersVLC);
            CEUseriVLC = zeros(1,UsersVLC);
            CSVLCaux = zeros(1,UsersVLC);
            for f=1:UsersVLC
                Waux = Wotimoi;
                Waux(:,f)=[];
                CLUseriVLC(f) = (log2((sum(alphaki.*((hVLCi(:,f))'*Wotimoi).^2)+1)/(sum(betaki.*(((hVLCi(:,f))'*Waux).^2))+ 1))/2);
                CEUseriVLC(f) = (log2((sum(alphae.*((gEveVLC)'*Wotimoi).^2)+1)/(sum(betae.*(((gEveVLC)'*Waux).^2))+ 1))/2);
                CSVLCaux(f)=max(CLUseriVLC(f)-CEUseriVLC(f),0);
            end
            sumVLCi = sum(CSVLCaux);
        else
            noiseVLC = 0;
            alphaki = (2^(2*log(1-(-1)))./(2*pi*exp(1).*noiseVLC));
            betaki = (1./(3*noiseVLC));
            Wotimoi=0;
            CSVLCaux = 0;
            sumVLCi = 0;
        end
        if  sumVLCi + sum(CSRFaux) > sumVLC + sumRF 
            alocUsers(i)=0;
            [CSUseri]=CalculoCS(alocUsers,hRFi,heRF,hVLCi,gEveVLC,Wotimoi,Votimoaux,noiseRF,betaki,betae,alphaki,alphae,nUsers);
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
        alphaki = (2^(2*log(1-(-1)))./(2*pi*exp(1).*noiseVLC));
        betaki = (1./(3*noiseVLC));
        [Wotimoi]=algoritmoVLC(noiseVLC,noiseEveVLC,rate,I_DC,gamma,eta,hVLCi,gEveVLC,UsersVLC,nLeds);
        CLUseriVLC=zeros(1,UsersVLC);
        CEUseriVLC=zeros(1,UsersVLC);
        CSVLCaux = zeros(1,UsersVLC);
            for z=1:UsersVLC
                W = Wotimoi;
                W(:,z)=[];
                CLUseriVLC(z) = (log2((sum(alphaki.*((hVLCi(:,z))'*Wotimoi).^2)+1)/(sum(betaki.*(((hVLCi(:,z))'*W).^2))+ 1))/2);
                CEUseriVLC(z) = (log2((sum(alphae.*((gEveVLC)'*Wotimoi).^2)+1)/(sum(betae.*(((gEveVLC)'*W).^2))+ 1))/2);
                CSVLCaux(z)=max(CLUseriVLC(z)-CEUseriVLC(z),0);
            end
        [Votimoaux]=algoritmoRFIterativo(noiseRF,rate,gammaSrf,hRFi,heRF,UsersRF,nAntenas);
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
        alphaki = (2^(2*log(1-(-1)))./(2*pi*exp(1).*noiseVLC));
        betaki = (1./(3*noiseVLC));
        [Wotimoi]=algoritmoVLC(noiseVLC,noiseEveVLC,rate,I_DC,gamma,eta,hVLCi,gEveVLC,UsersVLC,nLeds);
        CLUseriVLC=zeros(1,UsersVLC);
        CEUseriVLC=zeros(1,UsersVLC);
        CSVLCaux = zeros(1,UsersVLC);
            for z=1:UsersVLC
                W = Wotimoi;
                W(:,z)=[];
                CLUseriVLC(z) = (log2((sum(alphaki.*((hVLCi(:,z))'*Wotimoi).^2)+1)/(sum(betaki.*(((hVLCi(:,z))'*W).^2))+ 1))/2);
                CEUseriVLC(z) = (log2((sum(alphae.*((gEveVLC)'*Wotimoi).^2)+1)/(sum(betae.*(((gEveVLC)'*W).^2))+ 1))/2);
                CSVLCaux(z)=max(CLUseriVLC(z)-CEUseriVLC(z),0);
            end
        [Votimoaux]=algoritmoRFIterativo(noiseRF,rate,gammaSrf,hRFi,heRF,UsersRF,nAntenas);
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
    [Wotimoi]=algoritmoVLC(noiseVLC,noiseEveVLC,rate,I_DC,gamma,eta,hVLCi,gEveVLC,UsersVLC,nLeds);
    [Votimoaux]=algoritmoRFIterativo(noiseRF,rate,gammaSrf,hRFi,heRF,UsersRF,nAntenas);
    [CSUserfinal]=CalculoCS(alocUsers,hRFi,heRFtil,hVLCi,gEveVLCtil,Wotimoi,Votimoaux,noiseRF,betaki,betae,alphaki,alphae,nUsers);
    SumCS = sum(sumVLC)+sum(sumRF);
    sumcomp = sum(CSUserfinal);
else
    hUsers=hUsersi;
    CSRF=CSRFi;
    CSVLC=CSVLCi;
    CSUser=CSUseri;
end
end
sumcompara(S) = sumcomp;
CSalocs(S,:)=CSUserfinal;
Indalocs(S,:)=alocUsers;
end
IndTotal(((q-1)*Sa+1:q*Sa),:)=Indalocs;
CSTotal(((q-1)*Sa+1:q*Sa),:) = CSalocs; 
sumexp(q,:) = sumcompara;
end
%save sumAlg2users(q).mat sumTotal
%save IndAlg2users(q).mat IndTotal

%save sumAlg7users4.mat sumcomp
%save CSTotalAlg7users4.mat CSUseri
%save IndTotalAlg7users3.mat alocUsers

save IndAlgpot20VLC30.mat IndTotal
save CSTotalAlgpot20VLC30.mat CSTotal
save sumAlgpot20VLC30.mat sumexp
toc
