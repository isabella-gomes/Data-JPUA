function [Wzf]=ZF(nLeds,nUsers,hVLC,hEveVLC,I_DC)
if nUsers<nLeds
    alpha = 0;
else
    alpha=nUsers/I_DC;
end

B = inv(hVLC'*hVLC + alpha.*eye(nUsers))*hVLC';
Baux = hVLC*inv(hVLC'*hVLC + alpha.*eye(nUsers));
blinha = sqrt(I_DC/(trace(B*Baux)));
Waux =(blinha)*Baux;
if norm(Waux,inf)<=I_DC
    Wzf=Waux;
else
    Wzf = zeros(nLeds,nUsers);
end
% cvx_clear;
%     cvx_begin quiet
%     variable Wzf(nLeds,nUsers)
%     find Wzf;
%     subject to
%     hVLC'*Wzf == eye(nUsers);
%     norm(Wzf,inf)<=I_DC;
%     cvx_end
%     
%     if isnan(Wzf) | cvx_optval == -Inf
%         Wzf = zeros(nLeds,nUsers);
%     end
end