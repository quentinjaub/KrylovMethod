% Script flop par rapport a p :

%Test Block Gmres
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Chargement des matrices
[A,~] = Choix_A();
%[A,matname]=bfw398a;
%  load('mat1.mat','A');
% load('mat1.mat','b');

n = size(A,1);
p = input('Nombre de colonnes de B ?( > 10 ): ');

precond = input('Do you want a pre-conditioning ? (0/1) : ')
if (precond == 0)
  M1 = eye(size(A));
  M2 = eye(size(A));
else 
  [M1,M2] = Choix_M(A);
end

for i = 1:p
    b(:,i) = random('bino',n,0.5)*[1:n] + random('bino',n,0.5);
    %b = [b b+i];
end
x0= ones(size(A,1),p);
tol = 1e-13;
maxit=size(A,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fonction calculant le nombre de flops de chaque méthode
FlopGmres =@(i,p) nbFlop( i,A,0,p );
FlopBlock = @(i,p) nbFlop( i,A,1,p );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xmatlab=[];
resMatlab = [];
TimeBLock = [];
TimeMyGmres = [];
itBLock = [];
itMyGmres = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:p
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % My GMRES Block
    tic
    [ x,flag,relres,iter,resvec ] = gmresblock( A,b(:,1:i),tol,maxit, M1, M2,x0(:,1:i) );
    MyT = toc;
    TimeBLock(i) = MyT;
    itBLock(i) = iter*FlopBlock(iter,i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % My GMRES
    IterMyGmres = 0;
    tic
    for j = 1:i
        [Xreal,~,~,miter,resvecreal] = MyGMRES(A,b(:,j),x0(:,j),tol,maxit,M1,M2);
        Xmatlab = [Xmatlab , Xreal];
        resMatlab = [resMatlab , resvecreal ];
        IterMyGmres = IterMyGmres+miter;
    end

    MatlabT = toc;
    TimeMyGmres(i) = MatlabT;
    itMyGmres(i) = IterMyGmres*FlopGmres(IterMyGmres,i);
%     disp([' Temps d execution de My GMRES Block : ' num2str(MyT)])
%     disp([' Temps d execution de Matlab GMRES : ' num2str(MatlabT)])
%   disp([' Nombre d iterations de My GMRES : ' num2str(iter)])
%     disp([' Nombre d iterations de Matlab GMRES : ' num2str(IterMyGmres)])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Affichage
%     disp(['Norme de la difference entre le resultat de matlab et le notre : ' num2str(norm(Xmatlab-x,'fro'))]);
    Xmatlab = [];
end
figure(2)
plot([1:p],TimeBLock,'r');
hold on
plot([1:p],TimeMyGmres,'b');
legend('GMRES Block','My Gmres')
title('Time = f(p)');
figure(3)
plot([1:p],itBLock,'r');
hold on
plot([1:p],itMyGmres,'b');

legend('GMRES Block','My Gmres')
title('Iterations = f(p)')
