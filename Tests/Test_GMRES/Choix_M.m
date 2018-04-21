function [M1,M2] = Choix_M(A)
  choix = menu( 'Test Gradient Conjugue Préconditionné', ...
              'diagnonal','cholesky incomplet sans fill-in',...
              'LU incomplete ');
  switch choix
      case 1,
        % Diagonal
        M1 = diag(diag(A));
        M2 = eye(size(A));
      case 2,
        % Cholesky Incomplet sans remplissage
        M1 = ichol(A);
        M2 = M1';
      case 3,
         [M1,M2] = ilu( A );
    end
end