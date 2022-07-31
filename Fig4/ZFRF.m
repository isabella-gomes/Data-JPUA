function [Vzf]=ZFRF(nAntenas,nUsers,hbRF,heRF,gammaSrf)
noiseRF = 10^(-100/10);
rate = 1;
gammak = 2^rate-1;
[Vaux] = ValorInicial(noiseRF,gammaSrf,gammak,hbRF,nUsers,nAntenas);
for i=1:nUsers
    sumn(i) = (norm(Vaux(:,i),2)).^2;
end
potInicio = sum(sumn);
%potInicio = gammaSrf;
if nUsers<nAntenas
    alpha = 0;
else
    alpha=nUsers/potInicio;
end
B = inv(hbRF'*hbRF + alpha.*eye(nUsers))*hbRF';
Baux = hbRF*inv(hbRF'*hbRF + alpha.*eye(nUsers));
blinha = sqrt(potInicio/(trace(B*Baux)));
Va = (blinha)*Baux;
for i=1:nUsers
    sumna(i) = (norm(Va(:,i),2)).^2;
end
if sum(sumna)<gammaSrf
    Vzf = Va;
else
    Vzf = zeros(nAntenas,nUsers);
end
end