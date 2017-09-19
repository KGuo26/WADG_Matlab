% function [V Vr Vs V1 V2 V3 id] = bern_tri(N,r,s)
% 
% % use equivalence between W&B and equispaced nodes - get ordering
% [re se] = EquiNodes2D(N); [re se] = xytors(re,se);
% Ve = bern_tri_b(N,re,se);
% for i = 1:size(Ve,2)
%    [val iid] = max(Ve(:,i)); 
%    id(i) = iid;
% %    id(i) = i;
% end
% 
% [V Vr Vs V1 V2 V3] = bern_tri_b(N,r,s);
% V  = V(:,id);
% Vr = Vr(:,id); Vs = Vs(:,id);
% V1 = V1(:,id); V2 = V2(:,id); V3 = V3(:,id);

function [V Vr Vs VL1 VL2 VL3] = bern_tri(N,r,s)

% barycentric version
L1 = -(r+s)/2; L2 = (1+r)/2; L3 = (1+s)/2;
dL1r = -.5; dL2r = .5; dL3r = 0;
dL1s = -.5; dL2s = 0; dL3s = .5;

sk = 1;
% for i = 0:N
%     for j = 0:N-i
%         k = N-i-j;
for k = 0:N
    for j = 0:N-k
        i = N-j-k;
        C=factorial(N)/(factorial(i)*factorial(j)*factorial(k));
        V(:,sk) = C*(L1.^i).*(L2.^j).*(L3.^k);
        
        dL1 = C*i*(L1.^(i-1)).*(L2.^j).*(L3.^k);
        dL2 = C*j*(L1.^(i)).*(L2.^(j-1)).*(L3.^k);
        dL3 = C*k*(L1.^(i)).*(L2.^j).*(L3.^(k-1));
        if i==0
            dL1 = zeros(size(dL1));
        end
        if j==0
            dL2 = zeros(size(dL2));
        end
        if k == 0
            dL3 = zeros(size(dL3));
        end
        Vr(:,sk) = dL1.*dL1r + dL2.*dL2r + dL3.*dL3r;        
        Vs(:,sk) = dL1.*dL1s + dL2.*dL2s + dL3.*dL3s;
        
        VL1(:,sk) = dL1;
        VL2(:,sk) = dL2;
        VL3(:,sk) = dL3;
        sk = sk + 1;
    end
end

return


