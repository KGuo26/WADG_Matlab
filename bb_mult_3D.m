  % polynomial multiplication with Bernstein polynomials. 
% should be directly extendable to 3D

clear

N = 3;
N2 = 1;
[re se te] = EquiNodes3D(N+N2); [re se te] = xyztorst(re,se,te);
[rq sq  tq wq] = tet_cubature(3*(N+N2));
Vq = bern_basis_tet(N,rq,sq,tq);
Vq2 = bern_basis_tet(N+N2,rq,sq,tq);

% we may wish to rescale the representations of the weighting functionb
% basis to reduce numerical roundoff. 
VM = bern_basis_tet(N2,rq,sq,tq);  

M = (Vq2'*diag(wq)*Vq2);
Pq = M \ (Vq2'*diag(wq));

Lvals = {};
Lids = {};
Lsum = 0;
for i = 1:size(VM,2)    
    
    % polynomial multiplication by basis function B^M_i
    L = Pq*(diag(VM(:,i))*Vq);
    L(abs(L)<1e-8) = 0;    
    Lsum = Lsum + L;
    
    % use fact that L has one entry per column -> (L*f)_{Lids{i}} = (Lvals{i}.*f)
    % save values of L and row ids
    Lvals{i} = L(find(L));
    [ir ic] = find(L);
    Lids{i} = ir;
end

% test polynomial multiplication
f = (1:size(Vq,2))';
g = (1:size(VM,2))';
fg_exact = Pq*((VM*g).*(Vq*f));

N3 = N + N2;
NMp = (N3+1)*(N3+2)*(N3+3)/6;
fg = zeros(NMp,1);
for i = 1:size(VM,2)
    iids = Lids{i};
    fg(iids) = fg(iids) + g(i)*Lvals{i}.*f;
end

% should be small
norm(fg_exact-fg)