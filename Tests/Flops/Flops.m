clear all;
close all;
n = 573;
nnzA = 3829;
m=[1:100];
p = 3;
figure(1)

for i = 1:100
    FlopGmres = 2*i*(i+1)*n + 2*i*nnzA-i*n+2*sum(([1:i]+1).*([1:i].^2)-(1/3)*([1:i].^3));
    plot(i,FlopGmres,'r*')
    hold on;
    FlopBlock = i*p*n^2 - 2*(p^3)*i/3 + 2*sum(p*([1:i]+1).*((p*[1:i]).^2)-(1/3)*(p*[1:i]).^3) + 2*p*n*i.*(i+1)
    plot(i,FlopBlock,'b*')
end
legend('My Gmres Flops (r)','Gmres Block Flops (b)')
title('Flops = f(m) pour p=3, n=573, nnz(A)=3829')
