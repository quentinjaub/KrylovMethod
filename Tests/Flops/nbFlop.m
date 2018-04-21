function [ F ] = nbFlop( i,A,flag,p )
%Retourne le nombre de flops
% flag = 0 -> Mygmres
% flag = 1 -> Block Gmres
    n = size(A,1);
    if (flag == 0)
        F = 2*i*(i+1)*n + 2*i*nnz(A)-i*n+4*sum(([1:i]+1).*([1:i].^2)-(1/3)*([1:i].^3));
    elseif(flag == 1)
        F = i*p*n^2 - 2*(p^3)*i/3 + 4*sum(p*([1:i]+1).*((p*[1:i]).^2)-(1/3)*(p*[1:i]).^3) + 2*p*n*i.*(i+1);
    else
        error('Flag');
    end

end

