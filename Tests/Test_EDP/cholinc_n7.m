function [ L ] = cholinc_n7( A, droptol )
%CHOLINC_N7 simule l'ancienne fonction Matlab cholinc
%           plus robuste pour la factorisation de Cholesky incomplète
%           avec threshold (parce que pas vrai factorisation de Cholesky !!)

[L,U] = ilu(A,struct('type','ilutp','droptol',droptol,'thresh',0));
R = diag(sqrt(abs(diag(U))))\U;
L = R';
end

