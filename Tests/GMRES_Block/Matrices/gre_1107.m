function [A,matname] = gre_1107;
% Load the matrix gre-1107
%
A = rdcoord('gre_1107.mtx');
if nargout == 2
  matname = 'GRE1107';
end
%
%end
