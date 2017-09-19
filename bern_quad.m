function [V Vr Vs] = bern_quad(N,r,s)

r = (1+r)/2; % convert to unit
s = (1+s)/2; % convert to unit

V = zeros(length(r),(N+1)^2);
Vr = zeros(length(r),(N+1)^2);
Vs = zeros(length(r),(N+1)^2);
sk = 1;
for j = 0:N
    for i = 0:N    
        V(:,sk) = bern_1D(N,i,r).*bern_1D(N,j,s);
        Vr(:,sk) = .5*d_bern_1D(N,i,r).*bern_1D(N,j,s);
        Vs(:,sk) = .5*bern_1D(N,i,r).*d_bern_1D(N,j,s);
        sk = sk + 1;
    end
end

function bi = bern_1D(N,i,r)

bi = nchoosek(N,i)*(r.^i).*(1-r).^(N-i);

function dbi = d_bern_1D(N,i,r)

if (i==0)
    dbi = -N*(1-r).^(N-1);
elseif (i==N)
    dbi = N*(r.^(N-1));
else
    dbi = nchoosek(N,i)*r.^(i - 1).*(1 - r).^(N - i - 1).*(i - N*r);
end

