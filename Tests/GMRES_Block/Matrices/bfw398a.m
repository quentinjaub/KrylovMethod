function [A,matname] =  bfw398a;
%
A = rdcoord('bfw398a.mtx');
if nargout == 2
  matname = 'BFW398A';
end
%
%end
