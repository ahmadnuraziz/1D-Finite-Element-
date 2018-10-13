function y = int_hat1_f(xL,xR)

%This function evaluate \int_{xL}^{xR} f(x)*\phi(x) dx,
% where \phi(x) = (x-xL)/(xR-xL),using the Simpson rule.

xM = 0.5*(xL+xR);
y = (xR-xL)/6.0*(f(xL)*hat1(xL,xL,xR) + 4*f(xM)*hat1(xM,xL,xR)...
    + f(xR)*hat1(xR,xL,xR));
return