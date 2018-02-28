function Wave2D_manufactured

Globals2D

M=1;
N = 4;
K1D = 16;
FinalTime = 0.5;

cfun = @(x,y) 1 + 0.5*sin(pi*x).*sin(pi*y); % smooth velocity
%cfun = @(x,y) 1;

% manufactured solution
%ffun = @(x,y,t)  pi*(-1./(1+0.5*sin(pi*x).*sin(pi*y))+2).*sin(pi*x).*sin(pi*y).*sin(pi*t);
ffun = @(x,y,t)  pi*(-1./(cfun(x,y))+2).*sin(pi*x).*sin(pi*y).*sin(pi*t);
%ffun = @(x,y,t) pi*sin(pi*x).*sin(pi*y).*sin(pi*t);
pfun = @(x,y,t) sin(pi*x).*sin(pi*y).*cos(pi*t);
ufun = @(x,y,t) -cos(pi*x).*sin(pi*y).*sin(pi*t);
vfun = @(x,y,t) -sin(pi*x).*cos(pi*y).*sin(pi*t);

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

StartUp2D;

% used for plotting
[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

% used for quadrature
Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V; 
xq = Vq*x; yq = Vq*y;

%construct Pq, Pq is the projection matrix for Nodal Basis
Pq=V*V'*Vq'*diag(wq);

%construct the matrix C
Cq=cfun(xq,yq);

%construct projection for M
VM = Vandermonde2D(M,r,s);
PqM = VM*VM'*Vq'*diag(wq);
CM = PqM*Cq;
PNM = BB_projection2D(N,M);
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

% outer time step loop
tstep = 0;

while (time<FinalTime)
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        timelocal = time + rk4c(INTRK)*dt;
        f=ffun(x,y,timelocal);
        [rhsp, rhsu, rhsv] = acousticsRHS2D_manu(p,u,v,f);
       
        % initiate and increment Runge-Kutta residuals
        %apply invM*M_c^2
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
       
        % update fields
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        p = p+rk4b(INTRK)*resp;
        
        %p1 = pfun(x,y,timelocal);
        %error = norm(p-p1,'fro')
        
    end
    
    % Increment time
    time = time+dt;
    tstep = tstep+1;
    
end

p1 = pfun(x,y,FinalTime);
error = norm(p-p1,'fro')
