function [ A , b ] = Choix_A(  )

disp('Choisissez la matrice A de depart : ')
disp('0 -> mat0 ')
disp('1 -> mat1 ')
disp('2 -> mat2 ')
disp('3 -> mat3 ')
num = input(' ? : ')
if (num == 0)
    load('mat0.mat','A');
    load('mat0.mat','b');
elseif(num == 1)
%Chargement des matrices
    load('mat1.mat','A');
    load('mat1.mat','b');
elseif (num == 2)
    load('mat2.mat','A');
    load('mat2.mat','b');
else
%Chargement des matrices
    load('mat3.mat','A');
    load('mat3.mat','b');
end

end

