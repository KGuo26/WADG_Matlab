% firstly we generate the interpolation points on each element, which 
clear

N=4;
K1D=16;

cfun = @(x,y) sin(pi*x).*sin(pi*y);

[Nv,VX,VY,K,EToV] = unif_tri_mesh(K1D);

Nfp = N+1; Np = (N+1)*(N+2)/2; Nfaces=3; NODETOL = 1e-12;
[x,y] = Nodes2D(N); [r,s] = xytors(x,y);
V = Vandermonde2D(N,r,s); invV = inv(V);
MassMatrix = invV'*invV;
va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)';
x = 0.5*(-(r+s)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc));
y = 0.5*(-(r+s)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc));


C = cfun(x,y); 

C_modal = inv(V)*C;

tol = 10^-3;

m=zeros(512,1);
for k=1:512
 v =abs(C_modal(:,k));

 
if(v(2)>tol|v(6)>tol)
    m(k)=1;
end

if(v(3)>tol | v(7)>tol | v(10)>tol)
    m(k)=2;
end

if(v(4)>tol|v(8)>tol|v(11)>tol|v(13)>tol)
    m(k)=3;
end

if(v(5)>tol|v(9)>tol|v(12)>tol|v(14)>tol|v(15)>tol)
    m(k)=4;
end

end
