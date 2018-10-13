function y = int_hat2_f(xL,xR)

%This function evaluate \int_{xL}^{xR} f(x)*\phi(x) dx,
% where \phi(x) = (xR-x)/(xR-xL),using the Simpson rule.

xM = 0.5*(xL+xR);
y = (xR-xL)/6.0*(f(xL)*hat2(xL,xL,xR) + 4*f(xM)*hat2(xM,xL,xR)...
    + f(xR)*hat2(xR,xL,xR));
return