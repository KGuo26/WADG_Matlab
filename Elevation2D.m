function [ E] = Elevation2D(A,B,N)
%generate the elevation matrix from N to N-1


for i=1:(N)
    for j=1:(N+1)
        if sum(abs(A(j,:)-(B(i,:)+[1 0])))<10^-10
            k1=j;
        end
        if sum(abs(A(j,:)-(B(i,:)+[0 1])))<10^-10
            k2=j;
        end        

    end
    E(i,[k1 k2])=[(B(i,1)+1)/N (B(i,2)+1)/N];
end
E=E';
end

