function  [fg] =  BB_mult_2D(f,g,N,N2)

[re se] = EquiNodes2D(N+N2); [re se] = xytors(re,se);
[rq sq wq] = Cubature2D(3*(N+N2));
Vq = bern_basis_tri(N,rq,sq);
Vq2 = bern_basis_tri(N+N2,rq,sq);

VM = bern_basis_tri(N2,rq,sq); 

M = (Vq2'*diag(wq)*Vq2);
Pq = M \ (Vq2'*diag(wq));

Lvals = {};
Lids = {};
for i = 1:size(VM,2)
    
    % polynomial multiplication by basis function B^M_i
    L = Pq*(diag(VM(:,i))*Vq);
    L(abs(L)<1e-8) = 0;    
    
    % use fact that L has one entry per column -> (L*f)_{Lids{i}} = (Lvals{i}.*f)
    % save values of L and row ids
    Lvals{i} = L(find(L));
    [ir ic] = find(L);
    Lids{i} = ir;
end

N3 = N + N2;
NMp = (N3+1)*(N3+2)/2;
fg = zeros(NMp,1);
for i = 1:size(VM,2)
    iids = Lids{i};
    fg(iids) = fg(iids) + g(i)*Lvals{i}.*f;
end



end