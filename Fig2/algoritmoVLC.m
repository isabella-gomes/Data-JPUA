function [Wotimo]=algoritmoVLC(noiseVLC,noiseEveVLC,rate,I_DC,gamma,eta,gVLC,gEveVLC,nUsers,nLeds)
alphak = (2^(2*log(1-(-1)))./(2*pi*exp(1).*noiseVLC))';
alphae = (2^(2*log(1-(-1)))./(2*pi*exp(1).*noiseEveVLC))';
betak = (1./(3*noiseVLC))';
betae = (1./(3*noiseEveVLC))';
%gammak = 0.5;
convergence = 0;
Lmax = 10;
m=0;
%W = gVLC*inv(gVLC'*gVLC);
W = (-1 + (1+1)*rand(nLeds,nUsers));
while convergence==0 && m<=Lmax
    m = m+1;
    Wtil = W;
    if nUsers==1
        p3til=sum(alphae'.*(gEveVLC'*Wtil).^2);
        cvx_clear;
          cvx_begin quiet
          variable W(nLeds,nUsers)
          variables r1(1,nUsers)  r3(1,nUsers)  p1(1,nUsers)  p3(1,nUsers) 
          maximize sum(r1-r3)
          subject to
                 p1ux = alphak'.*(power(gVLC'*Wtil,2)+2*Wtil'*gVLC*gVLC'*(W-Wtil));
             p1<=sum(p1ux);
             p3>=sum(alphae'.*power(gEveVLC'*W,2));
             r1<=log(1+p1)/(2*log(2));
             r3>=log2(1+p3til)/2+(p3-p3til)/(2*log(2*(1+p3til)));
             r1>=rate;
          norm(W,inf)<=I_DC;
        cvx_end;
    else
        p2til = zeros(1,nUsers);
        p3til = zeros(1,nUsers);
        for i=1:nUsers
              Waux = Wtil;
              Waux(:,i)=[];
              p2til(i) = sum(betak*((gVLC(:,i)'*Waux).^2));
              p3til(i)=sum(alphae'.*(gEveVLC'*Wtil(:,i)).^2);
        end
        cvx_clear;
          cvx_begin quiet
          variable W(nLeds,nUsers)
          variables r1(1,nUsers) r2(1,nUsers) r3(1,nUsers) r4(1,nUsers) p1(1,nUsers) p2(1,nUsers) p3(1,nUsers) p4(1,nUsers)
          maximize sum(r1-r2-r3+r4)
          subject to
          p1ux = cvx(zeros(1,nUsers));
          p4aux = cvx(zeros(1,nUsers));
          for i=1:nUsers
              Wa = W;
              Wa(:,i)=[];
              Waux1 = Wtil;
              Waux1(:,i)=[];
              for j=1:(nUsers)
                 p1ux(j) = alphak'.*(power(gVLC(:,i)'*Wtil(:,j),2)+2*Wtil(:,j)'*gVLC(:,i)*gVLC(:,i)'*(W(:,j)-Wtil(:,j)));
              end
             for k=1:(nUsers-1)
                 p4aux(k) = betae'.*(power(gEveVLC'*Waux1(:,k),2)+2*Waux1(:,k)'*gEveVLC*gEveVLC'*(Wa(:,k)-Waux1(:,k)));
             end
             p1(i)<=sum(p1ux);
             p2(i)>=sum(betak'.*power(gVLC(:,i)'*Wa,2));
             p3(i)>=sum(alphae'.*power(gEveVLC'*W(:,i),2));
             p4(i)<=sum(p4aux);
             r1(i)<=log(1+p1(i))/(2*log(2));
             r2(i)>=log2(1+p2til(i))/2+(p2(i)-p2til(i))/(2*log(2*(1+p2til(i))));
             r3(i)>=log2(1+p3til(i))/2+(p3(i)-p3til(i))/(2*log(2*(1+p3til(i))));
             r4(i)<=log(1+p4(i))/(2*log(2));
             r1(i)-r2(i)>=rate;
          end
          norm(W,inf)<=I_DC;
        cvx_end;
    end
    if cvx_optval == -Inf | isnan(W)
        Wotimo = zeros(nLeds,nUsers);
        break
    end
        %disp(W)
        if norm(W-Wtil)/norm(W)<= 0.5
            convergence=1;
            Wotimo=W;
        elseif m>Lmax
            Wotimo=W;
        end
end
end