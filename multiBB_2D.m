function h=multiBB_2D(f,g,m,n)
%% f and g are the coefficient vectors of ordered Bernstein basis
%% m is the degree of f, n is the degree of g

%decide the polynomial with higher degree 
if m >= n
    f1 = f;
    g1 = g;
    m1 = m;
    n1 = n;
else
    f1 = g;
    g1 = f;
    m1 = n;
    n1 = m; 
    
end

l = m1 + n1;

h = zeros((l+1)*(l+2)/2,1);

for i=1:(n1+1)
    for j = 1:(n1+2-i)
        for k=1:(m1+1)
            for l=1:(m1+2-k)
               h(0.5*(i+k-2)*(2*n1+2*m1+5-i-k)+j+l-1)=h(0.5*(i+k-2)*(2*n1+2*m1+5-i-k)+j+l-1)+g1(0.5*(i-1)*(2*n1+4-i)+j)*f1(0.5*(k-1)*(2*m1+4-k)+l);
            end
        end
    end
end
    




end