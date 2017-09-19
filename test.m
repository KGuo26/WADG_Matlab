f1=[1 1 1 1 1 1 1 1 1 1]';
g1=[1 1 1 1]';
m1=2;
n1=1;


l = m1+n1;
h = zeros((l+1)*(l+2)*(l+3)/6,1);
c = 0;
c1 = 0;
c2 = 0;

for i1 = 1:(n1+1)
    for i2 = 1:(n1+2-i1)
        for i3 = 1:(n1+3-i1-i2)
            for j1 = 1:(m1+1)
                for j2 = 1:(m1+2-j1)
                    for j3 = 1:(m1+3-j1-j2)
                        
                        for k1 = 1:(i1+j1-2)
                            c = c + (l+3-k1)*(l+2-k1)/2;
                        end
                        c = c+(2*l+5-i2-j2)*(i2+j2-2)/2;
                        
                        
                        for k2 = 1:(i1-1)
                            c1 =c1+ (n1+3-k2)*(n1+2-k2)/2;
                        end
                        c1 = c1 + (2*n1+4-i2)*(i2-1)/2;
                
                        for k3 = 1:(j1-1)
                            c2 =c2+ (m1+3-k3)*(m1+2-k3)/2;
                        end
                        c2 = c2 + (2*m1+4-j2)*(j2-1)/2;
                        
                        h(c+i3+j3-1) = h(c+i3+j3-1) + g1(c1+i3)*f1(c2+j3);
                        
                        c=0;
                        c1=0;
                        c2=0;
                        
                    end
                end
            end
        end
    end
end
