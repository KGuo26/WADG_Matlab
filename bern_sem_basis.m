% wedge basis 
% function V = bern_sem_basis(N,r,s,t)
%     Vtri = bern_tri(N,r(:),s(:));
% %     [r2D s2D] = Nodes2D(N); [r2D s2D] = xytors(r2D,s2D);
% %     VWB = Vandermonde2D(N,r2D,s2D); Vtri = Vandermonde2D(N,r(:),s(:))/VWB;
% r1D = JacobiGL(0,0,N);
% VSEM = Vandermonde1D(N,r1D(:));
% V1D = Vandermonde1D(N,t(:))/VSEM; % sem in vertical direction
% 
% Np = (N+1)*(N+2)/2;
% 
% sk = 1;
% for i = 1:N+1
%     for j = 1:Np
%         V(:,sk) = Vtri(:,j).*V1D(:,i);
%         sk = sk + 1;
%     end
% end

function [V Vr Vs Vt V1 V2 V3] = bern_sem_basis(N,r,s,t)

% [Vtri Vrtri Vstri V1tri V2tri V3tri] = bern_tri(N,r(:),s(:));

[r2D s2D] = Nodes2D(N); [r2D s2D] = xytors(r2D,s2D);
VWB = Vandermonde2D(N,r2D,s2D);
Vtri = Vandermonde2D(N,r(:),s(:))/VWB;
[Vrtri Vstri] = GradVandermonde2D(N,r(:),s(:));
Vrtri = Vrtri/VWB; Vstri = Vstri/VWB;
V1tri = Vrtri*0; V2tri = Vrtri*0; V3tri = Vrtri*0; % dummy arrays

r1D = JacobiGL(0,0,N);
VSEM = Vandermonde1D(N,r1D(:));
V1D = Vandermonde1D(N,t(:))/VSEM; % sem in vertical direction
Vt1D = GradVandermonde1D(N,t(:))/VSEM; % sem in vertical direction

Nptri = (N+1)*(N+2)/2;

sk = 1;
for i = 1:N+1
    for j = 1:Nptri
        V(:,sk)  = Vtri(:,j).*V1D(:,i);
        Vr(:,sk) = Vrtri(:,j).*V1D(:,i);
        Vs(:,sk) = Vstri(:,j).*V1D(:,i);
        Vt(:,sk) = Vtri(:,j).*Vt1D(:,i);
        V1(:,sk) = V1tri(:,j).*V1D(:,i);
        V2(:,sk) = V2tri(:,j).*V1D(:,i);
        V3(:,sk) = V3tri(:,j).*V1D(:,i);
        sk = sk + 1;
    end
end