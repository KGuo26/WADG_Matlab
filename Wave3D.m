clear

Globals3D;

N = 4;
FinalTime = 1;
cfun = @(x,y,z) ones(size(x));
%cfun = @(x,y,z) 1 + 0.5*sin(pi*x).*sin(pi*y).*sin(pi*z); % smooth velocity
%cfun = @(x,y,z) (1 + .5*sin(2*pi*x).*sin(2*pi*y) + (y > 0)); % piecewise smooth velocity

%pfun = @(x,y,z,t) cos(pi*x).*cos(pi*y).*cos(pi*z).*cos(sqrt(3)*pi*t);

 pfun = @(x,y,z,t)   sin(pi*x).*sin(pi*y).*sin(pi*z).*cos(pi*t);
 ufun = @(x,y,z,t)  -cos(pi*x).*sin(pi*y).*sin(pi*z).*sin(pi*t);
 vfun = @(x,y,z,t)  -sin(pi*x).*cos(pi*y).*sin(pi*z).*sin(pi*t);
 wfun = @(x,y,z,t)  -sin(pi*x).*sin(pi*y).*cos(pi*z).*sin(pi*t);
 %ffun = @(x,y,z,t)  sqrt(3)*pi*(1-1./cfun(x,y,z)).*cos(pi*x).*cos(pi*y).*cos(pi*z).*sin(sqrt(3)*pi*t);
 ffun = @(x,y,z,t)   pi*2*sin(pi*x).*sin(pi*y).*sin(pi*z).*sin(pi*t);

%generate 3D mesh
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGmsh3D('cube2.msh');

StartUp3D;


Dr(abs(Dr)<1e-8) = 0; 

Nq = 2*N+1;

[rq sq tq wq] = tet_cubature(Nq); % integrate u*v*c

Vq = Vandermonde3D(N,rq,sq,tq)/V;

xq = Vq*x; yq = Vq*y; zq=Vq*z;
 
Pq=V*V'*Vq'*diag(wq);
 
Cq = cfun(xq,yq,zq);

Lift3D(N, r, s, t);
% %% Adptive m on each element
% VMq = Vandermonde3D(M,rq,sq,tq);
% CqM = VMq*VMq'*diag(wq)*Cq;



%% initial condition
p = pfun(x,y,z,0);
u = ufun(x,y,z,0);
v = vfun(x,y,z,0);
w = wfun(x,y,z,0);

% p = pfun(x,y,z,0);
% u = zeros(Np,K);
% v = zeros(Np,K);
% w = zeros(Np,K);



%%
time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resw = zeros(Np,K); resp = zeros(Np,K);

CN = (N+1)*(N+2)*(N+3)/6; % trace inequality constant
CNh = max(CN*max(Fscale(:)));
dt = 2/CNh;

%dt = 0.00390625; %coincide with the time step in occa
% outer time step loop
tstep = 0;




while (time<FinalTime)
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        timelocal = time + rk4c(INTRK)*dt;
        
        f=ffun(x,y,z,timelocal);
        
        [rhsp, rhsu, rhsv, rhsw] = acousticsRHS3D_WADG(p,u,v,w,f);
        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        resw = rk4a(INTRK)*resw + dt*rhsw;
                
        % update fields
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        w = w+rk4b(INTRK)*resw;
        p = p+rk4b(INTRK)*resp;
        
    
    end
    
    time = time+dt; tstep = tstep+1;

    
     
     
end

p_exact = pfun(x,y,z,FinalTime);

p_exact_quadrature = pfun(xq,yq,zq,FinalTime);

p_quadrature = Vq * p;

[d1,d2]=size(p_quadrature);

error_accumulation = 0;

for j1=1:d2
    for j2=1:d1
        err = p_quadrature(j2,j1)-p_exact_quadrature(j2,j1);
        error_accumulation = error_accumulation + err*err*wq(j2)*J(1,j1);
    end
end

error_l2 = sqrt(error_accumulation)
%error_fro = norm(p-p_exact,'fro'); 



function [rhsp, rhsu, rhsv, rhsw] = acousticsRHS3D_WADG(p,u,v,w,f)

Globals3D;

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);
dw = zeros(Nfp*Nfaces,K); dw(:) = w(vmapP)-w(vmapM);

% evaluate upwind fluxes
ndotdU = nx.*du + ny.*dv + nz.*dw;

% % Impose reflective boundary conditions (p+ = -p-)
% ndotdU(mapB) = 0;
% dp(mapB) = -2*p(vmapB);

% Impose Dirichlet BCs
ndotdU(mapB) = 0;
p(vmapB) = 0;
dp(mapB) = 0;


tau = 1;
fluxp =  tau*dp - ndotdU;
fluxu =  (tau*ndotdU - dp).*nx;
fluxv =  (tau*ndotdU - dp).*ny;
fluxw =  (tau*ndotdU - dp).*nz;

pr = Dr*p; ps = Ds*p; pt = Dt*p;

dpdx = rx.*pr + sx.*ps + tx.*pt;
dpdy = ry.*pr + sy.*ps + ty.*pt;
dpdz = rz.*pr + sz.*ps + tz.*pt;

divU = Dr*(u.*rx + v.*ry + w.*rz) + Ds*(u.*sx + v.*sy + w.*sz)+ Dt*(u.*tx + v.*ty + w.*tz);

% compute right hand sides of the PDE's
rhsp =  -divU + f + LIFT*(Fscale.*fluxp)/2.0;
rhsu =  -dpdx + LIFT*(Fscale.*fluxu)/2.0;
rhsv =  -dpdy + LIFT*(Fscale.*fluxv)/2.0;
rhsw =  -dpdz + LIFT*(Fscale.*fluxw)/2.0;

rhsp = Pq*(Cq.*(Vq*rhsp));

return;

end


