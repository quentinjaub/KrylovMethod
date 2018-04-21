function [A,matname] =  bwm200;
%
A = rdcoord('bwm200.mtx');
if nargout == 2
  matname = 'BWM200';
end
%
%end
