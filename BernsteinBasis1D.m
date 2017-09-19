function [ B ] = BernsteinBasis1D( N )
k=1;
for i=0:N
    B(k,:)=[N-i,i];
k=k+1;
end


end