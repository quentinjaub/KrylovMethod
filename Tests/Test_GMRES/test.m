clear all;
close all;

%Chargement de la matrice
[ A , b ] = Choix_A();

x0= ones(size(A,2),1);
tol = 1e-14;
maxit=size(A,1);
precond = input('Do you want a pre-conditioning ? (0/1) : ')
if (precond == 0)
  M1 = eye(size(A));
  M2 = eye(size(A));
else 
  [M1,M2] = Choix_M(A);
end

%My GMRES
[ x,flag,relres,iter,resvec ] = MyGMRES( A,b,x0,tol,maxit, M1, M2 );
%Matlab GMRES
[Xreal,~,~,~,resvecreal] = gmres(A,b,[],tol,maxit,M1,M2,x0);
norm(x-Xreal)
figure(2)
semilogy(resvec/norm(M2\(M1\b)))

hold on
semilogy( resvecreal/norm(M2\(M1\b) ),'r--')
legend('MyGMRES','GMRES Matlab')
semilogy([1 iter],[tol tol],'k--')


