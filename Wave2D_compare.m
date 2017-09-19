clear
kd=[4 8 16 32 64];
for i=1:5 

Globals2D
N=4;
%K1D = 16;
K1D = kd(i);
c_flag = 0;
FinalTime = 0.5;
%cfun = @(x,y) ones(size(x));
cfun = @(x,y) 1+0.5*sin(pi*x).*cos(pi*y); % smooth velocity
%cfun = @(x,y) (1 + .5*sin(2*pi*x).*sin(2*pi*y) + (y > 0)); % piecewise smooth velocity

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

StartUp2D;

%% for the plotting
[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y; % get quadrature points on each element

%% generate the quadrature point 
Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V; 
xq = Vq*x; yq = Vq*y;


%% construct Pq, the projection to degree N
Pq=V*V'*Vq'*diag(wq);

%{
%% construct Pq2, the projection to degree m
m=4;
[x4,y4] = Nodes2D(m); [r4,s4] = xytors(x4,y4);
Pq2 = Vandermonde2D(m,r4,s4) * Vandermonde2D(m,rq,sq)'* diag(wq); 
invV1 = inv(Vandermonde2D(m,r4,s4));
%}
%% construct projection matrices for degree 1
V1 = Vandermonde2D(1,rq,sq);
V2 = Vandermonde2D(1,r,s);
Pq1 = V1 * V1' *diag(wq); %Pq1 * Cq will give the projected function values at quadrature points

%% construct the matrix C, function values at quadrature points
Cq = cfun(xq,yq);
Cq1 = Pq1*Cq;
%{
%% porject the Cq to degree 4 and convert to modal basis 
C_modal = invV1 * Pq2 * Cq;


%% find degree m for each element
mtol = 10^-3;
m = zeros(K,1);

for k=1:K
    v = abs(C_modal(:,k)); 
    if(v(2)>mtol || v(6)>mtol)
        m(k)=1;
    end

    if(v(3)>mtol || v(7)>mtol || v(10)>mtol)
        m(k)=2;
    end

    if(v(4)>mtol|| v(8)>mtol || v(11)>mtol || v(13)>mtol)
        m(k)=3;
    end

    if(v(5)>mtol || v(9)>mtol || v(12)>mtol || v(14)>mtol || v(15)>mtol)
        m(k)=4;
    end
end
%disp(m)

%% using m, we project Cq1 to degree 1 if m=1; If m=0, we set the corresponding column of Cq to be a constant
Cq1 = zeros(size(Cq,1),K);
for k=1:K
    i = m(k);
    v = C_modal(:,k);
    if(i==0)
        Cq1(:,k) = v(1)*ones(size(Cq,1),1);
    end
    if(i==1)
        Cq1(:,k) = Pq1 * Cq(:,k);
    end
end
%}
%% initial condition

x0 = 0; y0 = .1;
p = exp(-25*((x-x0).^2 + (y-y0).^2));
u = zeros(Np, K);
v=zeros(Np,K);


%% construct the comparison
p1 = p;
u1 = u;
v1 = v;

time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);

resu1 = resu; resv1 = resv; resp1 = resp;

% compute time step size
CN = (N+1)*(N+2)/2; % trace inequality constant
CNh = max(CN*max(Fscale(:)));
dt = 2/CNh;

%% outer time step loop
tstep = 0;

while (time<FinalTime)
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        timelocal = time + rk4c(INTRK)*dt;
        
        [rhsp, rhsu, rhsv] = acousticsRHS2D_WADG(p,u,v);
        [rhsp1, rhsu1,rhsv1] = acousticsRHS2D_WADG1(p1,u1,v1);
        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        
        resp1 = rk4a(INTRK)*resp1 + dt*rhsp1;
        resu1 = rk4a(INTRK)*resu1 + dt*rhsu1;
        resv1 = rk4a(INTRK)*resv1 + dt*rhsv1;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        p = p+rk4b(INTRK)*resp;
        
        u1 = u1+rk4b(INTRK)*resu1;
        v1 = v1+rk4b(INTRK)*resv1;
        p1 = p1+rk4b(INTRK)*resp1;
    end
    %{
    if 1 && nargin==0 && mod(tstep,10)==0
        clf
        vv = Vp*p1;
        plot3(xp,yp,vv,'.');
        axis equal
        axis tight
        colorbar
        title(sprintf('time = %f',time))
        drawnow
    end
   
    if 1 && nargin==0 && mod(tstep,10)==0
       E = invV * (p-p1);
       ae = norm(E.*J,'fro');
       re = ae/norm(p.*J,'fro');
       disp(['  at time =',num2str(time),'   absolute difference=',num2str(ae),  '   relative difference=',num2str(re)]);
       
    
    end
    %}
    
    time = time+dt; tstep = tstep+1;
end
    E = invV * (p-p1);
    ae(i) = norm(E.*J,'fro');
    %re = ae/norm(p.*J,'fro');
    %disp(['  at time =',num2str(time),'   absolute difference=',num2str(ae),  '   relative difference=',num2str(re)]);
       
end 
l = [1:5];
loglog(kd,ae)
title('Plot of difference as mesh varies ')
xlabel('K1D')
ylabel('difference')