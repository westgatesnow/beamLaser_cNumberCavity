%This file is to solve the transcendental equation
%   kappa*a^2=phi*sin^2(g*a*tau/2) 
%as a fuction of tauList.

%Input: tauList, dim = (1, N)
%Output: idaList, dim = (1, N)

function aList = getAList(tauList, kappa, rabi, density)

aList = tauList.*0;
for i = 1:size(tauList,2)
    aList(i) = getAValue(tauList(i),kappa,rabi,density);
end
    
