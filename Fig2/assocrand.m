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
%aloc = randi([1 length(resultaux)],1,1);
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
    %CSVLC=0;
    CSRF=zeros(1,length(CSUser(alocUsers==0)));
    %CSRF=0;
    hUsers=zeros(nLeds,nUsers);
    hVLC=zeros(nLeds,length(CSUser(alocUsers==1)));
    %hVLC = zeros(nLeds,1);
    %hRF = zeros(nAntenas,1);
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
%noiseVLC = real(2.*gamma.*echarge.*eta.*hVLC'.*I_DC.*B + 4.*pi.*echarge.*Ar.*gamma.*Xamb.*(1-cos(Psic)).*B+i_amp^2.*B);
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
      CLUser(i) = log2((sum(alphak.*((hVLC(:,j))'*Wotimo).^2)+1)/(sum(betak.*(((hVLC(:,j))'*W).^2))+ 1))/2;
      CEUser(i) = log2((sum(alphae.*((gEveVLCtil)'*Wotimo).^2)+1)/(sum(betae.*(((gEveVLCtil)'*W).^2))+ 1))/2;
      CSVLC(j)=max(CLUser(i)-CEUser(i),0);
      j=1+j;
    elseif alocUsers(i)==0
        V = Votimo;
        V(:,k)=[];
        CLUser(i) = log2(1+(((hRF(:,k)'*Votimo(:,k)).^2))./(sum((((hRF(:,k)'*V)).^2))+ noiseRF));
        CEUser(i) = log2(1+(((heRFtil'*Votimo(:,k)).^2))./(sum(((heRFtil'*V).^2))+noiseRF));
        CSRF(k)=max(CLUser(i)-CEUser(i),0);
        k=1+k;
    end
      CSUser(i) = max(CLUser(i)-CEUser(i),0);
end
%sumcompara(p) = sumcomp;
IndTotal(q,:)=alocUsers;
CSTotal(q,:) = CSUser; 
sumexp(q,:) = sum(CSUser);
end
%save sumAlg2users(q).mat sumTotal
%save IndAlg2users(q).mat IndTotal

%save sumAlg7users4.mat sumcomp
%save CSTotalAlg7users4.mat CSUseri
%save IndTotalAlg7users3.mat alocUsers

save IndAlg5usersrand.mat IndTotal
save CSTotalAlg5usersrand.mat CSTotal
save sumAlg5usersrand.mat sumexp
toc
