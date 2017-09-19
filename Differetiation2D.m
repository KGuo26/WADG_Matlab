%Construct differential matrix D0 D1 D2

B=BernsteinBasis2D(4);

%construct D0
for j=1:15
    for i=1:15
        if sum(abs(B(j,:)-B(i,:)))<10^-10
            
            D0(i,j)=B(j,1);
            
        end
        if sum(abs(B(j,:)-B(i,:)-[1 -1 0]))<10^-10
           
            D0(i,j)=B(j,2)+1;
        
        end        

        if sum(abs(B(j,:)-B(i,:)-[1 0 -1]))<10^-10
            
            D0(i,j)=B(j,3)+1;
        
        end
    end
   
end

for j=1:15
    for i=1:15
        if sum(abs(B(j,:)-B(i,:)))<10^-10
            
             D1(i,j)=B(j,2);
        
        end
        
        if sum(abs(B(j,:)-B(i,:)-[-1 1 0]))<10^-10
            
            D1(i,j)=B(j,1)+1;
        
        end        

        if sum(abs(B(j,:)-B(i,:)-[0 1 -1]))<10^-10
            
            D1(i,j)=B(j,3)+1;
        
        end
    end
    
end

for j=1:15
    for i=1:15
        if sum(abs(B(j,:)-B(i,:)))<10^-10
            
            D2(i,j)=B(j,3);
        
        end
        
        if sum(abs(B(j,:)-B(i,:)-[-1 0 1]))<10^-10
            
            D2(i,j)=B(j,1)+1;
        
        end        

        if sum(abs(B(j,:)-B(i,:)-[0 -1 1]))<10^-10
       
            D2(i,j)=B(j,2)+1;
        
        end
    end
   
end