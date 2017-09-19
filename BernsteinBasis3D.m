function [ B ] = BernsteinBasis3D( N )

k=1;
for h = 0:N
    for l = 0:N-h
        for j=0:N-h-l
        i=N-h-l-j;
        B(k,:)=[i j l h];
        k=k+1;
    
        end
     end
end    
end

