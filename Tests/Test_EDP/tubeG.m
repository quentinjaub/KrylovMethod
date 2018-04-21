function [x,y]=be2G(bs,s)
%
nbs = 12;

if nargin==0,
  x=nbs; % nombre de segments frontière au total.
  return
end

% Segments droits
d=[
  0 0 0 0 0 0 % abscisse curviligne de départ
  1 1 1 1 1 1 % abscisse curviligne de fin
  1 1 1 1 1 1 % label du sous-domaine à gauche
  0 0 0 0 0 0 % label du sous-domaine à droite
];
%
% arcs de cercle extérieur
dc=[
  pi   pi   % abscisse curviligne de départ
  pi/2 pi/2 % abscisse curviligne de fin
  1    0    % label du sous-domaine à gauche
  0    1    % label du sous-domaine à droite
];
%
% sous-domaine 2
d2=[
  0   pi    % abscisse curviligne de départ
  pi  2*pi  % abscisse curviligne de fin
  2   2     % label du sous-domaine à gauche
  1   1     % label du sous-domaine à droite
];
% sous-domaine 3
d3=[
  0   pi    % abscisse curviligne de départ
  pi  2*pi  % abscisse curviligne de fin
  3   3     % label du sous-domaine à gauche
  1   1     % label du sous-domaine à droite
];
%
d = [d dc d2 d3];

% La ligne suivante force la mise sous forme vecteur ligne
% du vecteur bs en entrée
bs1=bs(:)';

if find(bs1<1 | bs1>nbs),
  error('Non existent boundary segment number')
end

if nargin==1,
  x=d(:,bs1);
  return
end

x=zeros(size(s));
y=zeros(size(s));
[m,n]=size(bs);
if m==1 & n==1,
  bs=bs*ones(size(s)); % dans le cas ou on n'a qu'un seul label de segment
                       % en entrée, c'est à dire que "bs" est un scalaire,
                       % on duplique ce scalaire autant de fois que de points
                       % attendus en sortie (nombre d'abscisses dans "s").
elseif m~=size(s,1) | n~=size(s,2),
  error('bs doit être un scalaire ou avoir la même taille que s');
end


% Enfin, calcul des coordonnées des points sur chacun des segments
% en question en fonction des abscisses curvilignes contenues dans
% le tableau "s":
%
% Demi-largeur extérieure:
  L = 1.0;
% Rayon de la zone de rafinement:
  R = 0.1;

if ~isempty(s),

% Bords extérieurs au sous-domaine 1 (segments de droite).
  js=1;  ii=find(bs==js);
  if length(ii)
     x(ii)=interp1([d(1,js),d(2,js)],[R L],s(ii));
     y(ii)=interp1([d(1,js),d(2,js)],[0 0],s(ii));
  end
  js=2;  ii=find(bs==js);
  if length(ii)
     x(ii)=interp1([d(1,js),d(2,js)],[L L],s(ii));
     y(ii)=interp1([d(1,js),d(2,js)],[0 L],s(ii));
  end
  js=3;  ii=find(bs==js);
  if length(ii)
     x(ii)=interp1([d(1,js),d(2,js)],[L  0],s(ii));
     y(ii)=interp1([d(1,js),d(2,js)],[L  L],s(ii));
  end
  js=4;  ii=find(bs==js);
  if length(ii)
     x(ii)=interp1([d(1,js),d(2,js)],[-L -L],s(ii));
     y(ii)=interp1([d(1,js),d(2,js)],[ 0 -L],s(ii));
  end
  js=5;  ii=find(bs==js);
  if length(ii)
     x(ii)=interp1([d(1,js),d(2,js)],[-L  0],s(ii));
     y(ii)=interp1([d(1,js),d(2,js)],[-L -L],s(ii));
  end
  js=6;  ii=find(bs==js);
  if length(ii)
     x(ii)=interp1([d(1,js),d(2,js)],[ 0  0],s(ii));
     y(ii)=interp1([d(1,js),d(2,js)],[-L -R],s(ii));
  end

%  Bords extérieurs au sous-domaine 1 (points anguleux).
  js=7; ii=find(bs==js);
  if length(ii)
     x(ii)=  R + R*cos(s(ii));
     y(ii)= -R + R*sin(s(ii));
  end

%  Bords arc de cercle exterieur
  js=8; ii=find(bs==js);
  if length(ii)
     x(ii)=  L*cos(s(ii));
     y(ii)=  L*sin(s(ii));
  end

%  Bords intérieurs sous-domaine 2.

  XC2 = -0.5;
  YC2 = -0.5;
  R2  =  0.3;

  js=9;  ii=find(bs==js);
  if length(ii)
     x(ii)= XC2 + R2*cos(s(ii));
     y(ii)= YC2 + R2*sin(s(ii));
  end
  js=10;  ii=find(bs==js);
  if length(ii)
     x(ii)= XC2 + R2*cos(s(ii));
     y(ii)= YC2 + R2*sin(s(ii));
  end

  XC2 =  0.5;
  YC2 =  0.5;
  R2  =  0.3;

  js=11;  ii=find(bs==js);
  if length(ii)
     x(ii)= XC2 + R2*cos(s(ii));
     y(ii)= YC2 + R2*sin(s(ii));
  end
  js=12;  ii=find(bs==js);
  if length(ii)
     x(ii)= XC2 + R2*cos(s(ii));
     y(ii)= YC2 + R2*sin(s(ii));
  end

end
