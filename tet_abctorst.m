function [r s t] = tet_abctorst(a,b,c)

% transform to tet
r = 0.5*(1+a).*.5.*(1-b).*(1-c)-1; 
s = 0.5*(1+b).*(1-c)-1; 
t = c;