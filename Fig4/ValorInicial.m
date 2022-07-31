function [Vtil] = ValorInicial(noiseRF,gammaSrf,gammak,hbob,nUsers,N)
noise = 3*noiseRF;
cvx_clear;
    cvx_begin quiet
    variable lambda(nUsers)
    maximize sum(lambda*(noise))
    subject to
    for j=1:nUsers
        H=zeros(N,N);
        for i=1:nUsers
        H = H + lambda(i)*hbob(:,i)*hbob(:,i)';
        end
        eye(N)+H -(1+1/gammak).*lambda(j).*hbob(:,j)*hbob(:,j)'== hermitian_semidefinite(N);
    end
    lambda>=0;
    cvx_end
      Vtil = zeros(N,nUsers);
     for j=1:nUsers
         for i=1:nUsers
         sommai = sum(lambda(i)*hbob(:,i)*hbob(:,i)');
         end
         aux = inv(eye(N)+ sommai)*hbob(:,j);
         Vtil(:,j)= sqrt(gammaSrf).*(aux/norm(aux));
     end
end