function [ B ] = BernsteinBasis2D( N )

k=1;
for l=0:N
 for i=N-l:-1:0
    j=N-i-l;
    B(k,:)=[i j l];
    k=k+1;
 end
end
end

