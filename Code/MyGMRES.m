%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithme GMRES avec preconditionnement
% Les parametres de cette fonctions seront :
    % 1. A, la matrice du systeme que l’on cherche a r ́esoudre
    % 2. b, le second membre de ce syst`eme
    % 3. x0, le vecteur initial
    % 4. tol, le seuil demandée
    % 5. maxit, le nombre maximum d’iterations
    % 6. M1, M2 : Matrices de pre conditionnement
%Les sorties sont : 
    % 1. x, la solution
    % 2. flag qui indique si la m ́ethode a converg ́e
    % 3. relres, la norme relative du rŕesidu (≡ η b N (x m ))
    % 4. iter, le nombre d’it ́erations
    % 5. resvec, le vecteur des normes des r ́esidus de chaque it ́eration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ x,flag,relres,iter,resvec ] = MyGMRES( A,b,x0,tol,maxit, M1, M2 )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialisation des variables
r0 = M2\(M1\(b-A*x0));
norm_r0 = norm(r0);
resvec = norm_r0;
beta = norm_r0;
norm_RHS = norm(M2\(M1\b));
relres = norm_r0/norm_RHS; %Residu/norme du second membre preconditionne
j=1;
V=r0/beta; % Vecteur Vj 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while (j<maxit && relres>tol)
    w=M2\(M1\(A*V(:,j)));
    %Arnoldi Orthogonalization
    for i=1:j
        h(i,j) = V(:,i)'*w; %Matrice h hessenberg
        w = w - h(i,j)*V(:,i);
    end
    h(j+1,j) = norm(w);
    %Augmentation de V
    V(:,j+1) = w/h(j+1,j);
    %e1 : premier vecteur de la base canonique
    e1=zeros(j+1,1);
    e1(1)=1;
    
    %Calcul Residu
    % QR factorization in order to avoid computation of x
    [Q,R] = qr(h); 
    g = beta*Q'*e1;
    %Calcul residu
    resvec(j+1) = abs(g(end));
    relres = resvec(j+1)/norm_RHS;
    j=j+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calcul de la solution finale x
yj = R(1:end-1,:)\g(1:end-1);
x=x0+V(:,1:j-1 )*yj;
relres = norm(b-A*x)/norm_RHS;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Flag de sortie
if(j==maxit)
    %Iterations max
    flag=0;
else
    %Convergence de l'algorithme
    flag=1;
end
iter=j;
end


