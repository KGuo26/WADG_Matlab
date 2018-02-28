clear
kd=[4 8 16 32 64];
h = 2./kd;
for i=1:length(kd)

Globals2D
N=7;
M=5;
%K1D = 16;
K1D = kd(i);
c_flag = 0;
FinalTime = 5.0;

cfun = @(x,y) 1 + 0.5*sin(pi*x).*sin(pi*y); % smooth velocity
% manufactured solution
ffun = @(x,y,t)  pi*(-1./(cfun(x,y))+2).*sin(pi*x).*sin(pi*y).*sin(pi*t);
pfun = @(x,y,t) sin(pi*x).*sin(pi*y).*cos(pi*t);
ufun = @(x,y,t) -cos(pi*x).*sin(pi*y).*sin(pi*t);
vfun = @(x,y,t) -sin(pi*x).*cos(pi*y).*sin(pi*t);

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

%% construct projection matrices for degree 1
V1 = Vandermonde2D(M,rq,sq);
V2 = Vandermonde2D(M,r,s);
Pq1 = V1 * V1' *diag(wq); %Pq1 * Cq will give the projected function values at quadrature points

%% construct the matrix C, function values at quadrature points
Cq = cfun(xq,yq);
Cq1 = Pq1*Cq;

%% initial condition
p = pfun(x,y,0);
u = ufun(x,y,0);
v = vfun(x,y,0);

time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);

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
        f=ffun(x,y,timelocal);
        [rhsp, rhsu, rhsv] = acousticsRHS2D_manu(p,u,v,f);
        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        p = p+rk4b(INTRK)*resp;
    end
    time = time+dt; tstep = tstep+1;
end
    p_final = pfun(x,y,FinalTime);
    ae(i) = norm(p-p_final);
end 
loglog(h,ae)
title('Plot of difference as mesh varies ')
xlabel('h')
ylabel('difference')
hold on 
%loglog(h,h.^3)
