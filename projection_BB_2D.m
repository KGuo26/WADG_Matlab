function [P] = projection_BB_2D(N,M)


Np = (N+1)*(N+2)/2;
NMp = (N+M+1)*(N+M+2)/2;

[r s] = Nodes2D(N+M); [r s] = xytors(r,s);
VMN = bern_basis_tri(N+M,r,s);
VN = bern_basis_tri(N,r,s);

TN = Vandermonde2D(N,r,s)\VN;
TMN = Vandermonde2D(N+M,r,s)\VMN;
EMN = VMN \ VN;

P = pinv(EMN); 




end