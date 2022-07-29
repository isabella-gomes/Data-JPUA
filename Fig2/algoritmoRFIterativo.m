function [Votimo]=algoritmoRFIterativo(noiseRF,rate,gammaSrf,hbRF,heRF,nUsers,nAntenas)
betak = (1./(noiseRF));
betae = (1./(noiseRF));
gammak = 2^rate-1;
convergence = 0;
Lmax = 10;
m=0;
[V] = ValorInicial(noiseRF,gammaSrf,gammak,hbRF,nUsers,nAntenas);
%V = hbRF*inv(hbRF'*hbRF);
%V = (-1+(1+1)*rand(nAntenas,nUsers));
while convergence==0 && m<=Lmax
    m = m+1;
    Vtil = V;
    if nUsers==1
        p3til=sum(betae.*(heRF'*V).^2);
        cvx_clear;
          cvx_begin quiet
          variable V(nAntenas,nUsers)
          variables r1(1,nUsers)  r3(1,nUsers)  p1(1,nUsers)  p3(1,nUsers) 
          maximize sum(r1-r3)
          subject to
             p1ux = betak'.*(power(hbRF'*Vtil,2)+2*Vtil'*hbRF*hbRF'*(V-Vtil));
             p1<=sum(p1ux);
             p3>=sum(betae'.*power(heRF'*V,2));
             r1<=log(1+p1)/(2*log(2));
             r3>=log2(1+p3til)/2+(p3-p3til)/(2*log(2*(1+p3til)));
             r1>=rate;
             sum(pow_pos(norm(V,2),2))<=gammaSrf;
        cvx_end;
    else
        for i=1:nUsers
              Vaux = Vtil;
              Vaux(:,i)=[];
              p2til(i) = sum(betak*((hbRF(:,i)'*Vaux).^2));
              p3til(i)=sum(betae.*(heRF'*V(:,i)).^2);
        end
        cvx_clear;
          cvx_begin quiet
          variable V(nAntenas,nUsers)
          variables r1(1,nUsers) r2(1,nUsers) r3(1,nUsers) r4(1,nUsers) p1(1,nUsers) p2(1,nUsers) p3(1,nUsers) p4(1,nUsers)
          maximize sum(r1-r2-r3+r4)
          subject to
          p1ux = cvx(zeros(1,nUsers));
          p4aux = cvx(zeros(1,nUsers));
          sumn=cvx(zeros(1,nUsers));
          for i=1:nUsers
              Va = V;
              Va(:,i)=[];
              Vaux1 = Vtil;
              Vaux1(:,i)=[];
              for j=1:(nUsers)
                 p1ux(j) = betak.*(power(hbRF(:,i)'*Vtil(:,j),2)+2*Vtil(:,j)'*hbRF(:,i)*hbRF(:,i)'*(V(:,j)-Vtil(:,j)));
              end
             for k=1:(nUsers-1)
                 p4aux(k) = betae.*(power(heRF'*Vaux1(:,k),2)+2*Vaux1(:,k)'*heRF*heRF'*(Va(:,k)-Vaux1(:,k)));
             end
             p1(i)<=sum(p1ux);
             p2(i)>=sum(betak.*power(hbRF(:,i)'*Va,2));
             p3(i)>=sum(betae.*power(heRF'*V(:,i),2));
             p4(i)<=sum(p4aux);
             r1(i)<=log(1+p1(i))/(log(2));
             r2(i)>=log2(1+p2til(i))+(p2(i)-p2til(i))/(log(2*(1+p2til(i))));
             r3(i)>=log2(1+p3til(i))+(p3(i)-p3til(i))/(log(2*(1+p3til(i))));
             r4(i)<=log(1+p4(i))/(log(2));
             r1(i)-r2(i)>=rate;
             sumn(i) = pow_pos(norm(V(:,i),2),2);
          end
          sum(sumn)<=gammaSrf;
          %norm(V,inf)<=gammaSrf;
        cvx_end;
        saux = zeros(1,nUsers);
        for i=1:nUsers
            saux(i) = norm(V(:,i),2).^2;
        end
        if sum(saux)>gammaSrf
            Votimo = ValorInicial(noiseRF,gammaSrf,gammak,hbRF,nUsers,nAntenas);
            break
        end
    end
    if isnan(V) | cvx_optval == -Inf
        V = ValorInicial(noiseRF,gammaSrf,gammak,hbRF,nUsers,nAntenas);
        %break
    end
        %disp(V)
        if norm(V-Vtil)/norm(V)<= 0.5
            convergence=1;
            Votimo=V;
        elseif m>Lmax
            Votimo=V;
        end
end
end