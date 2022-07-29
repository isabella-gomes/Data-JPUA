function [CSUser]=CalculoCS(alocUsers,hRFi,heRF,hVLCi,gEveVLC,Wotimoi,Votimoi,noiseRF,betaki,betae,alphaki,alphae,nUsers)
CLUser = zeros(1,nUsers);
CEUser = zeros(1,nUsers);
CSUser = zeros(1,nUsers);
j=1;
k=1;
for i=1:nUsers
    if alocUsers(i)==1
     W = Wotimoi;
     W(:,j)=[];
      CLUser(i) = abs(log2((sum(alphaki.*((hVLCi(:,j))'*Wotimoi).^2)+1)/(sum(betaki.*(((hVLCi(:,j))'*W).^2))+ 1))/2);
      CEUser(i) = abs(log2((sum(alphae.*((gEveVLC)'*Wotimoi).^2)+1)/(sum(betae.*(((gEveVLC)'*W).^2))+ 1))/2);
      j=1+j;
    elseif alocUsers(i)==0
        V = Votimoi;
        V(:,k)=[];
        CLUser(i) = log2(1+(((hRFi(:,k)'*Votimoi(:,k)).^2))./(sum((((hRFi(:,k)'*V)).^2))+ noiseRF));
        CEUser(i) = abs(log2(1+(((heRF'*Votimoi(:,k)).^2))./(sum(((heRF'*V).^2))+noiseRF)));
        k=1+k;
    end
end
CSUser = max(CLUser-CEUser,0);