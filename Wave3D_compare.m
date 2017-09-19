function Wave3D_compare

format long e

Globals3D

N = 3;
K1D = 16;
c_flag = 0;
FinalTime = 0.5;
%cfun = @(x,y,z) ones(size(x));
cfun = @(x,y) 1 + .5*sin(pi*x).*sin(pi*y); % smooth velocity
% cfun = @(x,y) (1 + .5*sin(2*pi*x).*sin(2*pi*y) + (y > 0)); % piecewise smooth velocity

%generate 3D mesh
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGmsh3D('cube1.msh');

StartUp3D;

Nq = 2*N+1;
[rq sq tq wq] = tet_cubature(Nq); % integrate u*v*c
Vq = Vandermonde3D(N,rq,sq,tq)/V;
xq = Vq*x; yq = Vq*y; zq=Vq*z;

Pq=V*V'*Vq'*diag(wq);

Cq = cfun(xq,yq);


%% Adptive m on each element
m = 3; % the largest degree that we can project

PMq = V*V'*Vq'*diag(wq);

c_nodal = cfun(xq,yq);

c_modal = inv(V) * PMq * c_nodal;

mD = zeros(K,1);

mtol=10^-2;

%using threshold to compute the degree
for i=1:K
    
    v =abs(c_modal(:,i)); 
  
    if(v(2)>mtol|| v(5)>mtol||v(11)>mtol)
        mD(i)=1;
    end

    if(v(3)>mtol || v(6)>mtol || v(8)>mtol|| v(12)>mtol|| v(14)>mtol|| v(17)>mtol)
        mD(i)=2;
    end

    if(v(4)>mtol|| v(7)>mtol|| v(9)>mtol|| v(10)>mtol|| v(13)>mtol|| v(15)>mtol|| v(16)>mtol|| v(18)>mtol|| v(19)>mtol|| v(20)>mtol)
        mD(i)=3;
    end
end



%disp(mD)
%% Generate the coefficient of c^2 based on BB basis
c_BB1 = inv(TB)*V*c_modal;
c_BB = zeros(Np,K);
for i =1:K
    
    k = mD(i);
    
    if(k==0 || k==1)
        kp = (k+1)*(k+2)*(k+3)/6;
    
        PMN = BB_projection3D(k,N-k);
     
        c_BB(1:kp,i)= PMN*c_BB1(:,i);
    else
        c_BB(:,i) = c_BB1(:,i);
    end
end

%% initial condition

p = cos(pi*x).*cos(pi*y).*cos(pi*z);
u1 = zeros(Np, K);
u2 = zeros(Np, K);
u3 = zeros(Np, K);

q = cos(pi*x).*cos(pi*y).*cos(pi*z);
v1 = zeros(Np, K);
v2 = zeros(Np, K);
v3 = zeros(Np, K);
%%
time = 0;

% Runge-Kutta residual storage
resu1 = zeros(Np,K); resu2 = zeros(Np,K); resu3 = zeros(Np,K); resp = zeros(Np,K);
resv1 = zeros(Np,K); resv2 = zeros(Np,K); resv3 = zeros(Np,K); resq = zeros(Np,K);
% compute time step size
CN = (N+1)*(N+2)*(N+3)/2; % trace inequality constant
CNh = max(CN*max(Fscale(:)));
dt = 2/CNh;
%dt = 0.00390625; %coincide with the time step in occa
% outer time step loop
tstep = 0;



while (time<FinalTime)
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        timelocal = time + rk4c(INTRK)*dt;
        
        [rhsp, rhsu1, rhsu2, rhsu3] = acousticsRHS3D_WADG(p,u1,u2,u3);
        [rhsq, rhsv1, rhsv2, rhsv3] = acousticsRHS3D_apdative(q,v1,v2,v3);
        % initiate and increment Runge-Kutta residuals
        resp  = rk4a(INTRK)*resp  + dt*rhsp;
        resu1 = rk4a(INTRK)*resu1 + dt*rhsu1;
        resu2 = rk4a(INTRK)*resu2 + dt*rhsu2;
        resu3 = rk4a(INTRK)*resu3 + dt*rhsu3;
        
        resq  = rk4a(INTRK)*resq  + dt*rhsq;
        resv1 = rk4a(INTRK)*resv1 + dt*rhsv1;
        resv2 = rk4a(INTRK)*resv2 + dt*rhsv2;
        resv3 = rk4a(INTRK)*resv3 + dt*rhsv3;
        % update fields
        u1 = u1+rk4b(INTRK)*resu1;
        u2 = u2+rk4b(INTRK)*resu2;
        u3 = u3+rk4b(INTRK)*resu3;
        p = p+rk4b(INTRK)*resp;
        
        v1 = v1+rk4b(INTRK)*resv1;
        v2 = v2+rk4b(INTRK)*resv2;
        v3 = v3+rk4b(INTRK)*resv3;
        q = q+rk4b(INTRK)*resq;
    end
    
    if 1 && nargin==0 && mod(tstep,50)==0
        E = inv(V)*(p-q);
        e = norm(E.*J,2);
        
        wJq = diag(wq)*(Vq*J);
        sum(sum(wJq.*(Vq*p - pex(xq,yq)).^2))
        fprintf("absolute error at time:%f is %f\n",time,e);
        fprintf("relative error at time:%f is %f\n",time,e/norm(p.*J,2));
    end
    
    time = time+dt; tstep = tstep+1;


    

end


end