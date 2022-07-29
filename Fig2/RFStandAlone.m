%Definir quantidade de usuários em VLC e RF inicialmente
%Definir a alocação inicial
clear all
tic
elements = {0:1, 0:1,0:1, 0:1,0:1, 0:1}; %cell array with N vectors to combine
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
aloc = perms([0 0 1 1]);
indice = randi(length(aloc));
alocUsers = aloc(indice,:);
%alocUsers = [0 0];
%alocUsers = randi([0,1],1,nUsers);
convergence=0;
iteracao =0;
UsersVLC=0;
UsersRF=0;
CSUser=zeros(1,nUsers);
    CLUser=zeros(1,nUsers);
    CEUser=zeros(1,nUsers);
    CSRF=zeros(1,length(CSUser(alocUsers==0)));
    %CSVLC=0;
    %CSRF=0;
    hUsers=zeros(nLeds,nUsers);
    %hVLC=zeros(nLeds,length(CSUser(alocUsers==1)));
    %hVLC = zeros(nLeds,1);
    hRF = zeros(nAntenas,1);
    k=1;
for i=1:nUsers
    if UsersRF<nAntenas && alocUsers(i)==0
        hUsers(:,i)=hbRF(:,i);
        hRF(:,k)=hbRF(:,i);
        UsersRF=UsersRF+1;
        k=k+1;
    end
end

if UsersRF ~= 0
[Votimo]=algoritmoRFIterativo(noiseRF,rate,gammaSrf,hRF,heRF,UsersRF,nAntenas);
end
j=1;
k=1;

for i=1:nUsers
    if alocUsers(i)==0
        V = Votimo;
        V(:,k)=[];
        CLUser(i) = log2(1+(((hRF(:,k)'*Votimo(:,k)).^2))./(sum((((hRF(:,k)'*V)).^2))+ noiseRF));
        CEUser(i) = log2(1+(((heRF'*Votimo(:,k)).^2))./(sum(((heRF'*V).^2))+noiseRF));
        CSRF(k)=max(CLUser(i)-CEUser(i),0);
        k=1+k;
    elseif alocUsers(i)==1
        CLUser(i) = 0;
        CEUser(i) = 0;
    end
      CSUser(i) = max(CLUser(i)-CEUser(i),0);
end
CSRFi = CSRF;
hUsersi = hUsers;
hRFi = hRF;
sumRF = sum(CSRF);
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
        if alocUsers(i)==1 && UsersRF==(nAntenas-2)
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
        %ind = aux;
        %ind = find(CSUseri==minVLC,1);
        for b = 1:length(hRFi)
            if hRFi(:,b)==hbRF(:,ind)
                hRFi(:,b)=hbRF(:,i);
                break
            end
        end
        [Votimoi]=algoritmoRFIterativo(noiseRF,rate,gammaSrf,hRFi,heRF,UsersRF,nAntenas);
        CLUseriRF=zeros(1,UsersRF);
        CEUseriRF=zeros(1,UsersRF);
        CSRFaux = zeros(1,UsersRF);
            for z=1:UsersRF
                V = Votimoi;
                V(:,z)=[];
                CLUseriRF(i) = log2(1+(((hRFi(:,k)'*Votimoi(:,k)).^2))./(sum((((hRFi(:,k)'*V)).^2))+ noiseRF));
                CEUseriRF(i) = log2(1+(((heRF'*Votimoi(:,k)).^2))./(sum(((heRF'*V).^2))+noiseRF));
                CSRFaux(z)=max(CLUseriRF(z)-CEUseriRF(z),0);
            end
        if sum(CSRFaux)>sumRF       
           alocUsers(i)=0;
           alocUsers(ind)=1;
           %[CSUseri]=CalculoCS(alocUsers,hRFi,heRF,hVLCi,gEveVLC,Wotimoi,Votimoaux,noiseRF,betaki,betae,alphaki,alphae,nUsers);
                %UsersRF=UsersRF-1;
           CSRFi = CSRFaux;
           sumRF=sum(CSRFi);
        else 
           for b = 1:length(hRFi)
            if hRFi(:,b)==hbRF(:,i)
                hRFi(:,b)=hbRF(:,ind);
                break
            end
           end
           CSUseri=CSaux;
        end
        end
end
if alocUsers==alocTest
    convergence=1;
    [Votimoi]=algoritmoRFIterativo(noiseRF,rate,gammaSrf,hRFi,heRF,UsersRF,nAntenas);
    CSRFfim = zeros(1,UsersRF);
    for z=1:UsersRF
                V = Votimoi;
                V(:,z)=[];
                CLUseriRF(i) = log2(1+(((hRFi(:,k)'*Votimoi(:,k)).^2))./(sum((((hRFi(:,k)'*V)).^2))+ noiseRF));
                CEUseriRF(i) = abs(log2(1+(((heRFtil'*Votimoi(:,k)).^2))./(sum(((heRFtil'*V).^2))+noiseRF)));
                CSRFfim(z)=max(CLUseriRF(z)-CEUseriRF(z),0);
    end
    sumRFfim = sum(CSRFfim);
    SumCS = sum(sumRFfim);
else
    %SumCS = 0;
    hUsers=hUsersi;
    CSRF=CSRFi;
end
end
%sumcompara(p) = sumcomp;
IndTotal(q,:)=alocUsers;
CSTotal(q,:) = CSRFfim; 
sumexp(q,:) = sum(CSRFfim);
end
%save sumAlg2users(q).mat sumTotal
%save IndAlg2users(q).mat IndTotal

%save sumAlg7users4.mat sumcomp
%save CSTotalAlg7users4.mat CSUseri
%save IndTotalAlg7users3.mat alocUsers

save IndAlg2usersRF1.mat IndTotal
save CSTotalAlg2RF1.mat CSTotal
save sumAlg2usersRF1.mat sumexp
toc
