function [gammaSrf,noiseRF,nUsers,nLeds,Ar,dv,Psi_k,phi,Psic,phi12,l,Ts,r,gamma,eta,echarge,Pot_Led,B,Xamb,i_amp,gVLC,gEveVLC,gEveVLCtil,noiseEveVLC,nAntenas,hbRF,heRF,heRFtil]=CSSNR
%clc
clear all  
Prf = 20;
gammaSrf = 10^(Prf/10)/1000;
noiseRF = (10^(-100/10))/1000; 
%S = length(gammaS);
%Pled = gammaS/nLeds;
%% Variáveis Canal VLC
nUsers = 2;
nLeds = 4;
Ar = 1e-4; %area PD
dv = 2.4; % Altura da sala
Psi_k = (2*pi)/3; %angulo de incidencia
phi = (2*pi)/3; %angulo de irradiancia
Psic = pi/3; % FOV
phi12 = pi/3; % angulo de meia intensidade do led
l = -log(2)/(log(cos(phi12))); %ordem de emissão lambertiana
Ts = 1; %ganho do filtro optico
r = 1.5; % indice de refração

%%Variáveis ruido VLC
gamma = 0.54; %responsitividade
eta = 0.44; %fator de conversão do LED
Pot_Led = 10^(40/10)/1000;
echarge = 1.6*10^(-19);
B = 10*10^6; %largura de banda
Xamb = 10.93; %fotocorrente luz ambiente
i_amp = 5*10^(-12); %preamplifier noise current density

Leds = [sqrt(2) sqrt(2);-sqrt(2) sqrt(2); -sqrt(2) -sqrt(2); sqrt(2) -sqrt(2)]; % posição dos Leds
Raio = dv.*tan(phi12);
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

%Usuários distribuidos uniformemente no espaço
    users = -2.5 + (2.5+2.5)*rand(nUsers,2); %Posicionando os usuários aleatóriamente
    Eve = -2.5 + (2.5+2.5)*rand(1,2);
    
    Deltarf = 0.1;
    Deltavlc = 1e-10; %variancia imperfect CSI Eve
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
    gVLC= real(State1.*(Ar./(distsv.^2)).*((l+1)/(2*pi)).*cos(phi)^l.*Ts.*(r^2/(sin(Psic)^2)).*cos(Psi_k));
    noiseVLC = 10^(-21);
    %noiseVLC = real(2.*gamma.*echarge.*eta.*gVLC'.*I_DC.*B + 4.*pi.*echarge.*Ar.*gamma.*Xamb.*(1-cos(Psic)).*B+i_amp^2.*B);
    [ine1,one1] = inpolygon(Eve(1),Eve(2),xLed1,yLed1); %canal do espião
    [ine2,one2] = inpolygon(Eve(1),Eve(2),xLed2,yLed2);
    [ine3,one3] = inpolygon(Eve(1),Eve(2),xLed3,yLed3);
    [ine4,one4] = inpolygon(Eve(1),Eve(2),xLed4,yLed4);
    eveVLC = repmat(Eve,nLeds,1);
    distsEve =  sqrt(((eveVLC(:,1)-Leds(:,1)).^2)+((eveVLC(:,2)-Leds(:,2)).^2));
    distEvevlc = sqrt(dv^2 + distsEve.^2);
    StateEve = [ine1;ine2;ine3;ine4];
    gEveVLC =real(StateEve.*(Ar./(distEvevlc.^2)).*((l+1)/(2*pi)).*cos(phi)^l.*Ts.*(r^2/(sin(Psic)^2)).*cos(Psi_k));
    gEveVLCtil = gEveVLC + sqrt(Deltavlc).*randn(length(gEveVLC),1);
    noiseEveVLC = 10^(-21);

    
    %% RF
RFfonte = [0,0];
Krice = 10.^(10/10);
phi = -1.8;
nAntenas = 4;



%% Canal RF - Bob e Eve
    distRF = sqrt(((users(:,1)-RFfonte(1)).^2)+((users(:,2)-RFfonte(2)).^2));
    distRFe = sqrt(((Eve(1)-RFfonte(1)).^2)+((Eve(2)-RFfonte(2)).^2));
    
    distRFv = sqrt(dv^2+distRF.^2);
    distRFve = sqrt(dv^2+distRFe.^2);
    OmegaRF = distRFv.^phi;
    OmegaRFe = distRFve.^phi;
    v = sqrt((Krice.*OmegaRF')./(1 + Krice)); %non-centrality parameter
    v = repmat(v,nAntenas,1);
    ve = sqrt((Krice.*OmegaRFe)./(1 + Krice));
    sigmab = sqrt(OmegaRF'./(2.*(1 + Krice))); %scale parameter
    sigmab = repmat(sigmab,nAntenas,1);
    sigmae = sqrt(OmegaRFe./(2.*(1 + Krice)));
    hbRF = sqrt(sigmab.^2.*(ncx2rnd(2,(v./sigmab).^2,nAntenas,nUsers)));
    heRF = sqrt(sigmae.^2.*(ncx2rnd(2,(ve./sigmae).^2,nAntenas,1)));
    heRFtil = heRF + sqrt(Deltarf/2)*(randn(length(heRF),1)+1i*randn(length(heRF),1));

end