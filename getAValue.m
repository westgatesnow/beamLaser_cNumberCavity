%This file is to solve the transcendental equation
%   kappa*a^2=phi*sin^2(g*a*tau/2) 
%as a fuction of tau.

function aValue = getAValue(tau, kappa, rabi, density)

myfun = @(a, t, k, r, d) d*r^2*t^2/4/k*sinc(r*a*t/2/pi)^2-1;
%matlab sinc(x) = sin(pi x)/(pi x)
k = kappa;
d = density;
r = rabi;
t = tau;
fun = @(a) myfun(a,t,k,r,d);
aValue = fzero(fun,2*pi/r/t);

    