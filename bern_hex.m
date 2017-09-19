function [V Vr Vs Vt] = bern_hex(N,r,s,t)

r = (1+r)/2; % convert to unit
s = (1+s)/2; % convert to unit
t = (1+t)/2; % convert to unit

Np = (N+1)^3;
V = zeros(length(r),Np);
Vr = zeros(length(r),Np);
Vs = zeros(length(r),Np);
Vt = zeros(length(r),Np);
sk = 1;
for k = 0:N
    for j = 0:N
        for i = 0:N
            
            V(:,sk) = bern_1D(N,i,r).*bern_1D(N,j,s).*bern_1D(N,k,t);
            Vr(:,sk) = .5*d_bern_1D(N,i,r).*bern_1D(N,j,s).*bern_1D(N,k,t);
            Vs(:,sk) = .5*bern_1D(N,i,r).*d_bern_1D(N,j,s).*bern_1D(N,k,t);
            Vt(:,sk) = .5*bern_1D(N,i,r).*bern_1D(N,j,s).*d_bern_1D(N,k,t);
            sk = sk + 1;
        end
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

