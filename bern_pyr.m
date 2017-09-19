function [V Vr Vs Vt Va Vb Vc] = bern_pyr(N,r,s,t)

a = 2*(r+1)./(1-t)-1;
b = 2*(s+1)./(1-t)-1;
c = t;
ids = abs(t-1)<1e-8;
a(ids) = -1;
b(ids) = -1;

dadr = 2./(1-t); 
dbds = 2./(1-t);
dadt = (1+a)./(1-t);
dbdt = (1+b)./(1-t);

sk = 1;
for k = 0:N
    for i = 0:N-k        
        for j = 0:N-k            
            V(:,sk) = bern(N-k,i,a).*bern(N-k,j,b).*bern(N,k,c);
            va = .5*d_bern(N-k,i,a).*bern(N-k,j,b).*bern(N,k,c);
            vb = .5*bern(N-k,i,a).*d_bern(N-k,j,b).*bern(N,k,c);
            vc = .5*bern(N-k,i,a).*bern(N-k,j,b).*d_bern(N,k,c);
            Vr(:,sk) = va.*dadr;
            Vs(:,sk) = vb.*dbds;
            Vt(:,sk) = va.*dadt + vb.*dbdt + vc;
            
            % for testing
            Va(:,sk) = va; Vb(:,sk) = vb; Vc(:,sk) = vc;            
            
            sk = sk + 1;            
        end
    end
end

function [bi] = bern(N,i,r)

r = (1+r)/2;
bi = nchoosek(N,i)*(r.^i).*(1-r).^(N-i);
% bi = (r.^i).*(1-r).^(N-i);

function dbi = d_bern(N,i,r)

if (i==0)
    dbi = -N*(1-r).^(N-1);
elseif (i==N)
    dbi = N*(r.^(N-1));
else
    dbi = nchoosek(N,i)*r.^(i - 1).*(1 - r).^(N - i - 1).*(i - N*r);
end