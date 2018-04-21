%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithme GMRES Par Bloc Avec Preconditionnement A Gauche
% Les parametres d’entree sont :
%   A, la matrice du syst`eme que l’on cherche `a resoudre,
%   B, le second membre de ce syst`eme,
%   tol, valeur seuil de l’erreur inverse `a utiliser pour detecter la convergence,
%   max it, le nombre d iterations maximum,
%   M1, M2, preconditionneur explicite en deux facteurs, `a utiliser pour accelerer la convergence
%   (signification identique `
%   a l’interface Matlab de gmres),
%   X0, la vecteur initial.
% Les resultats sont :
%   X, la solution calculee du syst`eme AX = B,
%   flag, statut de sortie de la routine ; ce param`etre sera  egal `a 0 en cas de convergence et `
%   a 1 en
%   cas de non convergence au bout de max it iterations.
%   relres, la norme relative du residu ( equivalent `a l erreur inverse η b N )
%   iter, le nombre d iterations,
%   resvec, le vecteur des normes des residus de chaque iteration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X, flag, relres, iter, resvec] = gmresblock(A, B, tol, maxit, M1, M2, X0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation
p = size(B,2);
R0 = M2\(M1\(B-A*X0));
norm_RHS = norm(M2\(M1\B),'fro');
[Q1,R] = qr(R0);
Vj = Q1(:,1:p);
R1 = R(1:p,1:p);
Stock_V = Vj;
resvec = norm(R0,'fro');
relres = resvec/norm_RHS;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j=0;
while (j<maxit && relres>tol)
    j = j+1;
    Wj = M2\(M1\(A*Vj));
    jp = (j-1)*p+1:j*p; %Vecteur pour inserer la nouvelle matrice Hij
    % Arnoldi
    for i = 1:j
        ip = (i-1)*p+1:i*p;
        H(ip,jp) = Stock_V(:,ip)'*Wj;
        Wj = Wj - Stock_V(:,ip)*H(ip,jp);
    end
	
    %Thin QR Factorization
    [Qw,Rw] = qr(Wj); 
    Vj = Qw(:,1:p);
    Stock_V = [Stock_V , Vj];
    H(j*p+1:(j+1)*p,jp) = Rw(1:p,1:p);
    
    % Base canonique
    E1=eye((j+1)*p,p);
    
    %Calcul residu
    [Qh,Rh] = qr(H); 
    g = Qh'*E1*R1;
    resvec(j+1) = norm(g(end,:),'fro');
    relres = resvec(j+1)/norm_RHS;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul de Y et X
yj = Rh(1:end-p,:)\g(1:end-p,:);
X=X0+Stock_V(:,1:end-p )*yj;

iter = j;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flag de sortie
if (j==maxit)
  flag = 0;
else
  flag = 1;
end

end
