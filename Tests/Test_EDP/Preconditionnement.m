% Résolution de l'équation elliptique
%   - nabla (c.grad u) + a.u = f

clear all
close all

% Les différents fichiers pour définir le problème d'EDP (voir explications)
% geometry
geom = 'tubeG';
% boundaries
boundary = 'tubeB';
% source
source = 'tubeF';

% commenter ces 2 lignes si pénurie de licences Matlab
figure(1);
pdegplot(geom), axis equal;

% Choix du niveau de raffinage : on effectuera autant de résolutions que de
% niveaux de raffinage
nR = input('Niveau de raffinage 0 <= nR < 4 : '); 

while (nR >= 4) || (nR < 0)
  nR = input('Niveau de raffinage 0 <= nR < 4 : ');
end

% Choix du préconditionneur pour le premier test
choix = menu( 'Test Gradient Conjugue Préconditionné', ...
              'sans','diagnonal','cholesky incomplet sans fill-in',...
              'cholesky incomplet avec tolérance','fin');

while choix < 5

  close all;
  
  % choix du seuil dans le cas de la factorisation incomplète de Cholesky avec
  % treshold
  if choix == 4
    DropTol = input('Drop Tolerance (réel positif) :  ') ;
  end

  % Création du maillage avec l'aide du fichier geom
  % commenter la ligne suivante si pénurie de licences Matlab
  [p,e,t] = initmesh(geom);

  % Boucle sur les niveaux de raffinage 
  for k = 0:nR

    % commenter cette section si pénurie de licences Matlab
    % -- début section

    % Raffinage
    if k > 0 
      [p,e,t] = refinemesh(geom, p, e, t);
    end

    % Dessin du maillage
    figure(2);
    pdemesh(p, e, t), axis equal
    xlabel(['number of triangles = ' num2str(size(t, 2))]);
    disp('fin construction du maillage : taper une touche');
    pause

    % Construction de la matrice de rigidité ainsi que du
    % second membre
    % problème résolu : - nabla(c.grad u ) + a.u = f
    % avec
    
    a = 0.0;
    
    % avec c variable suivant les sous-domaines Omega 1, 2 et 3
    
    c = setupC(p, t);
    
    % avec f donné par le fichier source
    % avec conditions aux limites données par le fichier boundary 
    %   Gamma 1 et Gamma 3 : Neumann homogènes
    %   Gamma 2 : Dirichlet u = 10
    %   Gamma 4 : Dirichlet u = 100
    
    [A,b]= assempde(boundary, p, e, t, c, a, source);

    % -- fin section

    % décommenter la section suivante si pénurie de licences Matlab
    % -- début section
    %switch k
    %  case 0
    %    load mat0;
    %  case 1
    %    load mat1;
    %  case 2
    %    load mat2;
    %  case 3
    %    load mat3;
    %  otherwise
    %    disp('impossible');
    %end
    % -- fin section

    % dimension du problème
    n = size(A, 1);

    % Définition des paramètres gouvernant l'arrêt de la méthode itérative
    tol = 1.e-10; 
    maxit = floor(n/2);

    % Construction du préconditionneur M1*M2 = M

    % tic
    timep = cputime;
    if (k==0)
        npcg = input('pcg or gmres ? (0/1) : ');
    end
    switch choix

      case 1,
        % Sans Préconditionnement
        M1 = eye(size(A));
        M2 = eye(size(A));
      case 2,
        % Diagonal
        M1 = diag(diag(A));
        M2 = eye(size(A));
      case 3,
        % Cholesky Incomplet sans remplissage
        M1 = ichol(A);
        M2 = M1';
      case 4,
        % Cholesky Incomplet avec treshold
         if (npcg == 0)
             M1 = cholinc_n7( A, DropTol );
             M2 = M1';
         else
             [M1,M2] = ilu( A );
         end

    end

    % tac
    timep = cputime - timep;

    if choix >= 3

      % Afficher la structure de la matrice dans la partie triangulaire
      % supérieure et celle du préconditionneur dans la partie triangulaire
      % inférieure.

      %Preconditionneur
      figure(7)
      spy(tril(M1))
      title('Structure Preconditionneur')
      % Matrix
      figure(8)
      spy(triu(A))
      title('Structure matrice A')

    end

    % Résolution du système préconditionné avec CG ou GMRES
    timer = cputime;
    
    % affectations pour permettre les affichages : à supprimer
    x = zeros(n, 1);
    iter = 100;
    % fin affectations à supprimer
    
    % À COMPLÉTER par l'appel à la fonction de résolution
    if (npcg == 0)
        %Conjuguate Gradient Method
        [x,flag,relres,iter,resvec] = pcg(A,b,tol,maxit,M1,M2);
    else
        %GMRES Solver
        [x,flag,relres,iter,resvec] = MyGMRES(A,b,x,tol,maxit, M1, M2);
    end
    timer= cputime - timer;
    
    % Affichage d'informations
    fprintf(' ------------------------------------------ \n');
    fprintf(' niveau de Raffinage : %5d \n', k);
    fprintf(' Taille du probleme : %5d \n', n);
    fprintf(' - Nb iterations : %4d \n' , iter);
    fprintf(' - CPU time pour la construction du préconditionneur : %3.1e s \n', timep);
    fprintf(' - CPU time pour la résolution : %3.1e s \n', timer);
    fprintf(' - CPU time : %3.1e s \n',timer + timep);

    % Dessin de la solution sur la géométrie
    % commenter les 4 lignes suivantes si pénurie de licences Matlab
    figure(4)
    Titre = [ 'Solution' ];
    pdeplot(p, e, t, 'xydata', x, 'title', Titre, 'colormap', 'jet', ...
            'mesh', 'off', 'contour', 'off', 'levels', 20), axis equal;
    fprintf(' ------------------------------------------ \n');
    disp('fin dessin solution : taper une touche');
    pause;

    % Afficher l'historique de convergence en les superposant pour les
    % différentes finesses de maillage

    % pour tracer les courbes de décroissance de la norme du résidu 
    % des différents maillages avec des couleurs différentes
    couleur = ['g', 'r', 'c', 'm'];

    figure(5)
    % À COMPLÉTER par le tracé de l'historique de convergence
    
    semilogy(resvec/norm(b),couleur(k+1))
    hold on
    title('Evolution of n_b(x_k)')
    
    disp('fin résolution pour ce maillage : taper une touche');
    pause

  end % for k, niveau de raffinage

  disp('nouveau calcul (avec autre préconditionneur ?)');
 
  % Choix du préconditionneur pour le test suivant
  choix = menu( 'Test Gradient Conjugue Préconditionné', ...
               'sans','diagnonal','cholesky incomplet sans fill-in',...
               'cholesky incomplet avec tolérance','fin');

end

close all;
