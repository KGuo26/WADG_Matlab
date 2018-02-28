format long
x = [0.776202 -1.714347 0.900284 -1.247643 -0.464629 ...
     1.231183 -2.158620 0.050615  0.306521 -0.402293 ...
     -0.993916 1.940297 1.289644 0.050972  0.919905 ...
     -0.114709 0.163441 -2.044624 -0.513900  0.542938]';
 
 rho=[1 1 1 1]';
 c=[1.6;-0.8;0.1333;-0.0048]; 
 y = mult_BB_3D(x, rho, 3,1);
 
 
 N = 3;
 M = 1;  
Np = (N+1)*(N+2)*(N+3)/6;
NMp = (N+M+1)*(N+M+2)*(N+M+3)/6;

[r s t]= Nodes3D(6); [r s t] = xyztorst(r,s,t);
VMN = bern_basis_tet(N+M,r,s,t);
VN = bern_basis_tet(N,r,s,t);

TN = Vandermonde3D(N,r,s,t)\VN;
TMN = Vandermonde3D(N+M,r,s,t)\VMN;
EMN = VMN \ VN;

PMN = pinv(EMN);

z = PMN*y;

for i = 0:N
    ENM{i+1} = bern_basis_tet(N+M,r,s,t)\bern_basis_tet(N-i,r,s,t);    
    EN{i+1} = bern_basis_tet(N,r,s,t)\bern_basis_tet(N-i,r,s,t);        
end

for i=1:4
    E{i} = bern_basis_tet(i,r,s,t)\bern_basis_tet(i-1,r,s,t); 
end

a = E{3}'*E{4}'*y;
b = E{2}'*a;
d = (c(3)*eye(4)+c(4)*E{1}*E{1}')*b;