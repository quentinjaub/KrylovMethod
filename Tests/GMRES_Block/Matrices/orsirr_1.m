function [A,matname] = orsirr_1;
A = rdcoord('orsirr_1.mtx');
if nargout == 2
  matname = 'ORSIRR_1';
end

