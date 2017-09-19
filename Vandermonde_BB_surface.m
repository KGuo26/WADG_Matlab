function [T0,T1,T2]=Vandermonde_BB_surface(r,s)

B=BernsteinBasis2D(4);
lambda0=-(r+s)/2;
lambda1=(1+r)/2;
lambda2=(1+s)/2;
lambda=[lambda0 lambda1 lambda2];

O1=[1 2 3 4 5];
O2=[5 9 12 14 15];
O3=[1 6 10 13 15];

for i=1:5
    for j=1:5
        T0(i,j)=factorial(4)/(factorial(B(O1(j),1))*factorial(B(O1(j),2))*factorial(B(O1(j),3)))*lambda0(O1(i))^B(O1(j),1)*lambda1(O1(i))^B(O1(j),2)*lambda2(O1(i))^B(O1(j),3);
    end
end

for i=1:5
    for j=1:5
        T1(i,j)=factorial(4)/(factorial(B(O2(j),1))*factorial(B(O2(j),2))*factorial(B(O2(j),3)))*lambda0(O2(i))^B(O2(j),1)*lambda1(O2(i))^B(O2(j),2)*lambda2(O2(i))^B(O2(j),3);
    end
end

for i=1:5
    for j=1:5
        T2(i,j)=factorial(4)/(factorial(B(O3(j),1))*factorial(B(O3(j),2))*factorial(B(O3(j),3)))*lambda0(O3(i))^B(O3(j),1)*lambda1(O3(i))^B(O3(j),2)*lambda2(O3(i))^B(O3(j),3);
    end
end

end


