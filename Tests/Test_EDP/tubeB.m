function [q,g,h,r]=be2G(p,e,u,time)
%
% Exécutable MATLAB spécifiant les conditions aux limites pour le problème
%   d'EDP
%
% Descripteur de conditions aux limites de Neumann homogènes :
b1 = [   ... 
 1       ... % dimension N of the system
 0       ... % number M of Dirichlet boundary conditions
 1       ... % length for the string representing q
 1       ... % length for the string representing g
 '0'     ... % scalar function q
 '0'     ... % scalar function g
]'*1;  l1 = length(b1);

% Descripteur de conditions aux limites de Neumann généralisées :
b2 = [   ... 
 1       ... % dimension N of the system
 0       ... % number M of Dirichlet boundary conditions
 3       ... % length for the string representing q
 3       ... % length for the string representing g
 '0.0'     ... % scalar function q
 '-10'     ... % scalar function g
]'*1;  l2 = length(b2);

% Descripteur de conditions aux limites de Dirichlet homogènes :
b3 = [   ... 
 1       ... % dimension N of the system
 1       ... % number M of Dirichlet boundary conditions
 1       ... % length for the string representing q
 1       ... % length for the string representing g
 1       ... % length for the string representing h
 3       ... % length for the string representing r
 '0'     ... % scalar function q
 '0'     ... % scalar function g
 '1'     ... % scalar function h
 '100'     ... % scalar function r
]'*1;  l3 = length(b3);

% Descripteur de conditions aux limites de Dirichlet homogènes :
b4 = [   ... 
 1       ... % dimension N of the system
 1       ... % number M of Dirichlet boundary conditions
 1       ... % length for the string representing q
 1       ... % length for the string representing g
 1       ... % length for the string representing h
 3       ... % length for the string representing r
 '0'     ... % scalar function q
 '0'     ... % scalar function g
 '1'     ... % scalar function h
 '10'     ... % scalar function r
]'*1;  l4 = length(b4);


lb = max([ l1 l2 l3 l4 ]);  cl = blanks(lb)'*ones(1,8);
cl(1:l1, 1) = b1;  
cl(1:l4, 2) = b4; 
cl(1:l1, 3) = b1;
cl(1:l1, 4) = b1;
cl(1:l3, 5) = b3;  
cl(1:l1, 6) = b1; 
cl(1:l1, 7) = b1;
cl(1:l1, 8) = b1;


if any(size(u))
  [q,g,h,r]=pdeexpd(p,e,u,time,cl);
else
  [q,g,h,r]=pdeexpd(p,e,time,cl);
end

