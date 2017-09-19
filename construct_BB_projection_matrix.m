clear

% Polynomial degree reduction in the L2-norm equals best Euclidean
% approximation of B?zier coefficients.
% D. Lutterkort, J. Peters, U. Reif 1999

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

% degree elevation matrices 
L = zeros(Np);
for i = 0:N
    ENM{i+1} = bern_basis_tet(N+M,r,s,t)\bern_basis_tet(N-i,r,s,t);    
    EN{i+1} = bern_basis_tet(N,r,s,t)\bern_basis_tet(N-i,r,s,t);        
end

for i=1:3
    E{i} = bern_basis_tet(i,r,s,t)\bern_basis_tet(i-1,r,s,t); 
end


% 2D - the formula for the eigenvalue should work for 3D if you set d = 3
d = 3;
lami = @(N,i) 2*factorial(N)^2*factorial(d)./(factorial(N+i+d).*factorial(N-i))';
eig_ratio = @(N,M,i) lami(N-i,0:N-i)./lami(N+M,0:N-i);

% code to test how to construct degree elevation operators
if 0
    % make i,j indices for degree N and degree N+M Bernstein bases
    sk = 1;
    for j = 0:N
        for i = 0:N-j
            idiN(sk) = i;
            idjN(sk) = j;
            sk = sk + 1;
        end
    end
    sk = 1;
    for j = 0:N+M
        for i = 0:N+M-j
            idiNM(sk) = i;
            idjNM(sk) = j;
            sk = sk + 1;
        end
    end
    % choose index to test
    i = 2;
    ratio = eig_ratio(N,M,i);
    idsc = (idiNM+idjNM <= (N-i)).*(idiNM+idjNM+1);
    idsr = (idiN+idjN <= (N-i)).*(idiN+idjN+1);
    col_ids = find(idsc);
    ratio_ids = idsc(col_ids);
    row_ids = find(idsr);
    EE = zeros(Np,NMp);
    for ii = 1:length(row_ids)
        EE(row_ids(ii),col_ids(ii)) = ratio(ratio_ids(ii));
    end
    norm(EE-TN*(EN{i+1}*ENM{i+1}'/TMN),'fro')
end

% use fact that inv(T)*E*T is diagonal to conclude that exists cj s.t.
% pinv(EMN) = sum_{j=0}^N E^N_{N-j}*E^(N+M)_{N-j} u 
L = zeros(N+1);
for ii = 0:N
    lam = eig_ratio(N,M,ii);
    L(1:length(lam),ii+1) = lam;    
end
b = ones(N+1,1); 
c = L\b;

A = c(1)*eye(20)+E{3}*(c(2)*eye(10)+E{2}*(c(3)*eye(4)+c(4)*E{1}*E{1}')*E{2}')*E{3}';

% construct P using coefficients
P = zeros(Np,NMp);
for ii = 0:N
    P = P + c(ii+1)*EN{ii+1}*ENM{ii+1}';
end
norm(P-PMN,'fro') % check if they're the same


