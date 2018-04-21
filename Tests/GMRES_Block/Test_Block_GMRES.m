%Test Block Gmres
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Chargement des matrices
[A,~] = Choix_A();
n = size(A,1);
p = input('Nombre de colonnes de B ? : ');

for i = 1:p
    b(:,i) = random('bino',n,0.5)*[1:n] + random('bino',n,0.5);
end
x0= ones(size(A,1),p);
tol = 1e-13;
maxit=size(A,1);

precond = input('Do you want a pre-conditioning ? (0/1) : ')
if (precond == 0)
  M1 = eye(size(A));
  M2 = eye(size(A));
else 
  [M1,M2] = Choix_M(A);
end

Xmatlab=[];
resMatlab = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% My GMRES Block
tic
[ x,flag,relres,iter,resvec ] = gmresblock( A,b,tol,maxit, M1, M2,x0 );
MyT = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab GMRES
matlabIter = 0;
tic
for j = 1:size(b,2)
    [Xreal,~,~,miter,resvecreal] = gmres(A,b(:,j),[],tol,maxit,M1,M2,x0(:,j));
    Xmatlab = [Xmatlab , Xreal];
    resMatlab = [resMatlab , resvecreal ];
    matlabIter = matlabIter+miter;
end
MatlabT = toc;
disp([' Temps d execution de My GMRES Block : ' num2str(MyT)])
disp([' Temps d execution de Matlab GMRES : ' num2str(MatlabT)])
disp([' Nombre d iterations de My GMRES Block : ' num2str(iter)])
disp([' Nombre d iterations de Matlab GMRES : ' num2str(matlabIter)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Affichage
disp(['Norme de la difference entre le resultat de matlab et le notre : ' num2str(norm(Xmatlab-x,'fro'))]);
figure(2)
semilogy(resvec/norm(M2\(M1\b)))

hold on
legend('GMRES Block')
semilogy([1 iter],[tol tol],'k--')


