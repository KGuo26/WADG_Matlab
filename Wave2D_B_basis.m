function Wave2D


Globals2D

N = 4;
K1D = 16;
c_flag = 0;
FinalTime = 0.5;
cfun = @(x,y) ones(size(x));
%cfun = @(x,y) 1 + .5*sin(pi*x).*sin(pi*y); % smooth velocity
%cfun = @(x,y) (1 + .5*sin(2*pi*x).*sin(2*pi*y) + (y > 0)); % piecewise smooth velocity

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

StartUp2D;
%construct LIFT matrix
LIFT_Bernstein;
%construct D0 D1 D2
Differetiation2D;

[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);

Vp = Vandermonde2D(N,rp,sp)/V;

xp = Vp*x; yp = Vp*y;


%Nq = 2*N+1;
%[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
%Vq = Vandermonde2D(N,rq,sq)/V; 
%xq = Vq*x; yq = Vq*y;



%construct Pq
%Pq=V*V'*Vq'*diag(wq);

%construct the matrix C
%C=cfun(xq,yq);




%% initial condition

x0 = 0; y0 = .1;
p = exp(-25*((x-x0).^2 + (y-y0).^2));
u = zeros(Np, K);
v=zeros(Np,K);
%%

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
        
        [rhsp, rhsu, rhsv] = acousticsRHS2D_BBprojection(p,u,v,timelocal);
        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        p = p+rk4b(INTRK)*resp;
        
    end;
    
    if 1 && nargin==0 && mod(tstep,10)==0
        clf
        vv = Vp*p;
        plot3(xp,yp,vv,'.');
        axis equal
        axis tight
        colorbar
        title(sprintf('time = %f',time))
        drawnow
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
end



