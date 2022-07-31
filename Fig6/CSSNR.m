
function [Leds,gammaSrf,noiseRF,nUsers,nLeds,Ar,dv,Psi_k,phi,Psic,l,Ts,r,gamma,eta,echarge,Pot_Led,B,Xamb,i_amp,hbRF,heRF,heRFtil,nAntenas,Deltarf,Deltavlc]=CSSNR(phi12,users,Eve)

%clc
%clear all
PsRf = 20;          
gammaSrf = 10.^(PsRf/10)/1000;
noiseRF = (10^(-100/10)/1000); 
%S = length(gammaS);
%Pled = gammaS/nLeds;
%% Variáveis Canal VLC
nUsers = 5;
nLeds = 4;
Ar = 1e-4; %area PD
dv = 2.4; % Altura da sala
Psi_k = (2*pi)/3; %angulo de incidencia
phi = (2*pi)/3; %angulo de irradiancia
Psic = pi/3; % FOV
%phi12 = pi/6; % angulo de meia intensidade do led
l = -log(2)./(log(cos(phi12))); %ordem de emissão lambertiana
Ts = 1; %ganho do filtro optico
r = 1.5; % indice de refração
%users = -2.5 + (2.5+2.5)*rand(3,2); %Posicionando os usuários aleatóriamente
%Eve = -2.5 + (2.5+2.5)*rand(1,2);
%%Variáveis ruido VLC
gamma = 0.54; %responsitividade
eta = 0.44; %fator de conversão do LED
Pot_Led = 10^(30/10)/1000;
echarge = 1.610^(-19);
B = 10*10^6; %largura de banda
Xamb = 10.93; %fotocorrente luz ambiente
i_amp = 5*10^(-12); %preamplifier noise current density

Leds = [sqrt(2) sqrt(2);-sqrt(2) sqrt(2); -sqrt(2) -sqrt(2); sqrt(2) -sqrt(2)]; % posição dos Leds
% Raio = dv.*tan(phi12);
% th = 0:pi/50:2*pi;
% %Posicionando Cobertura de cada fonte
% xLed1 = Raio.*cos(th)+Leds(1,1);
% yLed1 = Raio.*sin(th)+Leds(1,2);
% 
% xLed2 = Raio.*cos(th)+Leds(2,1);
% yLed2 = Raio.*sin(th)+Leds(2,2);
% 
% xLed3 = Raio.*cos(th)+Leds(3,1);
% yLed3 = Raio.*sin(th)+Leds(3,2);
% 
% xLed4 = Raio.*cos(th)+Leds(4,1);
% yLed4 = Raio.*sin(th)+Leds(4,2);

%Usuários distribuidos uniformemente no espaço
    
 Deltarf = 0.1;
 Deltavlc = 1e-10; %variancia imperfect CSI Eve   
    
    %% RF
RFfonte = [0,0];
Kdb = 10;
Krice = 10.^(Kdb/10);
percurso = -1.8;
nAntenas = 4;

%% Canal RF - Bob e Eve
    distRF = sqrt(((users(:,1)-RFfonte(1)).^2)+((users(:,2)-RFfonte(2)).^2));
    distRFe = sqrt(((Eve(1)-RFfonte(1)).^2)+((Eve(2)-RFfonte(2)).^2));
    distRFv = sqrt(dv^2+distRF.^2);
    distRFve = sqrt(dv^2+distRFe.^2);
    OmegaRF = distRFv.^percurso;
    OmegaRFe = distRFve.^percurso;
    v = sqrt((Krice.*OmegaRF')./(1 + Krice)); %non-centrality parameter
    v = repmat(v,nAntenas,1);
    ve = sqrt((Krice.*OmegaRFe)./(1 + Krice));
    sigmab = sqrt(OmegaRF'./(2.*(1 + Krice))); %scale parameter
    sigmab = repmat(sigmab,nAntenas,1);
    sigmae = sqrt(OmegaRFe./(2.*(1 + Krice)));
    hbRF = sqrt(sigmab.^2.*(ncx2rnd(2,(v./sigmab).^2,nAntenas,nUsers)));
    heRF = sqrt(sigmae.^2.*(ncx2rnd(2,(ve./sigmae).^2,nAntenas,1)));
    heRFtil = heRF + sqrt(Deltarf/2)*(randn(length(heRF),1)+1i*randn(length(heRF),1));
    %rate=1;
    %[Votimo]=algoritmoRFIterativo(noiseRF,rate,gammaSrf,hbRF,heRF,nUsers,nAntenas)
end