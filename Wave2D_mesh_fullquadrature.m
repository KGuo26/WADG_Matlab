clear
Globals2D

kd=[4 8 16 32];
h=2./kd;
N = 4;
M = 1;

for i = 1:length(kd)

K1D = kd(i);
FinalTime = 1.0;
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
StartUp2D;

%% Set up wavespeed function
%cfun = @(x,y) ones(size(x));
cfun = @(x,y) 1+ sin(pi*x).*sin(pi*y); % smooth velocity
%cfun = @(x,y) (1 + .5*sin(2*pi*x).*sin(2*pi*y) + (y > 0)); % piecewise smooth velocity

%% generate quadrature points
Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
xq = Vq*x; yq = Vq*y;

%% construct the projection matrix for nodal basis Pq
Pq = V*V'*Vq'*diag(wq);

%% construct matrix Cq
Cq = cfun(xq,yq);

%% construct the projection matrix for cfun to degree M
VMq = Vandermonde2D(M,rq,sq);
CqM = VMq*VMq'*diag(wq)*Cq;

%% initial condition
x0 = 0; y0 = .1;
p = exp(-25*((x-x0).^2 + (y-y0).^2));
u = zeros(Np, K);
v=zeros(Np,K);

p_M = p;
u_M = u;
v_M = v;

time = 0;

%% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);
resu_M = resu; resv_M = resv; resp_M = resp;

%% compute time step size
CN = (N+1)*(N+2)/2; % trace inequality constant
CNh = max(CN*max(Fscale(:)));
dt = 2/CNh;

%% outer time step loop
tstep = 0;

while (time<FinalTime)
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        timelocal = time + rk4c(INTRK)*dt;
        
        [rhsp, rhsu, rhsv] = acousticsRHS2D_fullWADG(p,u,v);
        [rhsp_M, rhsu_M, rhsv_M] = acousticsRHS2D_adaptiveWADG(p_M,u_M,v_M);
        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        
        resp_M = rk4a(INTRK)*resp_M + dt*rhsp_M;
        resu_M = rk4a(INTRK)*resu_M + dt*rhsu_M;
        resv_M = rk4a(INTRK)*resv_M + dt*rhsv_M;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        p = p+rk4b(INTRK)*resp;
        
        u_M = u_M + rk4b(INTRK)*resu_M;
        v_M = v_M + rk4b(INTRK)*resv_M;
        p_M = p_M + rk4b(INTRK)*resp_M;
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
end

p_quadrature = Vq * p;
p_M_quadrature = Vq * p_M;

[d1,d2]=size(p_quadrature);
error_accumulation = 0;
for j1=1:d2
    for j2=1:d1
        err = p_quadrature(j2,j1)-p_M_quadrature(j2,j1);
        error_accumulation = error_accumulation + err*err*wq(j2)*J(1,j1);
    end
end

error_l2(i) = sqrt(error_accumulation);
error_fro(i) = norm(p-p_M,'fro'); 
end

function [rhsp, rhsu, rhsv] = acousticsRHS2D_adaptiveWADG(p,u,v)

Globals2D;

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);

% evaluate upwind fluxes
ndotdU = nx.*du + ny.*dv;

% Impose reflective boundary conditions (p+ = -p-)
ndotdU(mapB) = 0;
dp(mapB) = -2*p(vmapB);

tau = 1;
fluxp =  tau*dp - ndotdU;
fluxu =  (tau*ndotdU - dp).*nx;
fluxv =  (tau*ndotdU - dp).*ny;

pr = Dr*p; ps = Ds*p;
dpdx = rx.*pr + sx.*ps;
dpdy = ry.*pr + sy.*ps;
divU = Dr*(u.*rx + v.*ry) + Ds*(u.*sx + v.*sy);

% compute right hand sides of the PDE's
rhsp =  -divU + LIFT*(Fscale.*fluxp)/2.0;
rhsu =  -dpdx + LIFT*(Fscale.*fluxu)/2.0;
rhsv =  -dpdy + LIFT*(Fscale.*fluxv)/2.0;

rhsp = Pq*(CqM.*(Vq*rhsp));
return;
end

function [rhsp, rhsu, rhsv] = acousticsRHS2D_fullWADG(p,u,v)

Globals2D;

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);

% evaluate upwind fluxes
ndotdU = nx.*du + ny.*dv;

% Impose reflective boundary conditions (p+ = -p-)
ndotdU(mapB) = 0;
dp(mapB) = -2*p(vmapB);

tau = 1;
fluxp =  tau*dp - ndotdU;
fluxu =  (tau*ndotdU - dp).*nx;
fluxv =  (tau*ndotdU - dp).*ny;

pr = Dr*p; ps = Ds*p;
dpdx = rx.*pr + sx.*ps;
dpdy = ry.*pr + sy.*ps;
divU = Dr*(u.*rx + v.*ry) + Ds*(u.*sx + v.*sy);

% compute right hand sides of the PDE's
rhsp =  -divU + LIFT*(Fscale.*fluxp)/2.0;
rhsu =  -dpdx + LIFT*(Fscale.*fluxu)/2.0;
rhsv =  -dpdy + LIFT*(Fscale.*fluxv)/2.0;

rhsp = Pq*(Cq.*(Vq*rhsp));
return;
end


