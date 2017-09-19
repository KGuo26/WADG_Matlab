function [V Vr Vs Vt] = bern_wedge(N,r,s,t)

[V1D Vt1D] = bern_basis_1D(N,t(:));

[Vtri Vrtri Vstri] = bern_basis_tri(N,r(:),s(:));

Np = (N+1)*(N+2)/2;
V = zeros(length(r(:)),Np*(N+1));
Vr = zeros(length(r(:)),Np*(N+1));
Vs = zeros(length(r(:)),Np*(N+1));
Vt = zeros(length(r(:)),Np*(N+1));
sk = 1;
for j = 1:N+1    
    for i = 1:Np    
        V(:,sk) = V1D(:,j).*Vtri(:,i);
        Vr(:,sk) = V1D(:,j).*Vrtri(:,i);
        Vs(:,sk) = V1D(:,j).*Vstri(:,i);
        Vt(:,sk) = Vt1D(:,j).*Vtri(:,i);
        
        sk = sk + 1;
    end
end