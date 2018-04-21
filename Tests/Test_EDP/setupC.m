function a=setupC(p,t)

 nt = size(t,2);

 a = ones(2,nt);

 ii = pdesdt(t, 2);
 a(1,ii) = 1e-3;
 a(2,ii)=1.0;

 ii = pdesdt(t, 3);
 a(1,ii) = 1e+3;
 a(2,ii) = 1.0;

