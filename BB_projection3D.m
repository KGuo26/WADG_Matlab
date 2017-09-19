function [PMN] = BB_projection3D(N,M)

Np = (N+1)*(N+2)*(N+3)/6;
NMp = (N+M+1)*(N+M+2)*(N+M+3)/6;

[r s t]= Nodes3D(N+M); [r s t] = xyztorst(r,s,t);
VMN = bern_basis_tet(N+M,r,s,t);
VN = bern_basis_tet(N,r,s,t);

TN = Vandermonde3D(N,r,s,t)\VN;
TMN = Vandermonde3D(N+M,r,s,t)\VMN;
EMN = VMN \ VN;

PMN = pinv(EMN); 


end