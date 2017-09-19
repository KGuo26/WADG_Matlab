
l=1;
for i=0:3
    for j=0:3-i
        for k=0:3-i-j
            A(l,:)=[i,j,k];
            l=l+1;
        end
    end
end