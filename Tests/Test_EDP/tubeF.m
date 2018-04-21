function f=be2F(p,t,u,time)

 nt = size(t,2);

 f = zeros(1,nt);

 ii = pdesdt(t, 2);
 f(ii) = 0;

