function [T]=Vandermonde_Bernstein(N,r,s)

B=BernsteinBasis2D(N);
lambda0=-(r+s)/2;
lambda1=(1+r)/2;
lambda2=(1+s)/2;

for i=1:15
    for j=1:15
        T(i,j)=factorial(N)/(factorial(B(j,1))*factorial(B(j,2))*factorial(B(j,3)))*(lambda0(i)^B(j,1))*(lambda1(i)^B(j,2))*(lambda2(i)^B(j,3));
    end
end

end

