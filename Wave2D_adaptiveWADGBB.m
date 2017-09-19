function Wave2D


Globals2D

N = 4;
K1D = 16;
c_flag = 0;
FinalTime = 0.5;
%cfun = @(x,y) ones(size(x));
cfun = @(x,y) 1+0.5*sin(pi*x).*sin(pi*y); % smooth velocity
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


Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V; 
xq = Vq*x; yq = Vq*y;


%construct Pq
Pq=V*V'*Vq'*diag(wq);

%construct the matrix C
C = cfun(xq,yq);

C_m = cfun(x,y);

C_modal = invV * Pq*C;%C_m;

mtol = 10^-2;

m = zeros(K,1);

for k=1:K
  v =abs(C_modal(:,k)); 
  if(v(2)>mtol|v(6)>mtol)
    m(k)=1;
  end

  if(v(3)>mtol | v(7)>mtol | v(10)>mtol)
    m(k)=2;
  end

  if(v(4)>mtol|v(8)>mtol|v(11)>mtol|v(13)>mtol)
    m(k)=3;
  end

  if(v(5)>mtol|v(9)>mtol|v(12)>mtol|v(14)>mtol|v(15)>mtol)
    m(k)=4;
  end
end

C_BB = inv(T) * C_m;

Proj{1} = projection_BB_2D(0,4);
Proj{2} = projection_BB_2D(1,3);
Proj{3} = projection_BB_2D(2,2);
Proj{4} = projection_BB_2D(3,1);
Proj{5} = projection_BB_2D(4,0);

C_BB_proj = zeros(Np,K);

for k=1:K
    mk = (m(k)+1)*(m(k)+2)/2;
    C_BB_proj(1:mk,k) = Proj{m(k)+1}*C_BB(:,k);
end


%% initial condition

x0 = 0; y0 = .1;
p = exp(-25*((x-x0).^2 + (y-y0).^2));
u = zeros(Np, K);
v=zeros(Np,K);
%%

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

% outer time step loop
tstep = 0;

while (time<FinalTime)
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        timelocal = time + rk4c(INTRK)*dt;
        
        [rhsp, rhsu, rhsv] = acousticsRHS2D_WADG(p,u,v);
        
        [rhsp1, rhsu1, rhsv1] = acousticsRHS2D_adaptiveBBprojection(p1,u1,v1,C_BB_proj,m);
        
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
    
    
    vv  = invV * p;
    vv1 = invV * p1;
    
    dv = abs(vv-vv1);
    
    error = 0;
    for i=1:K
        error = error + 0.039*sum(dv(:,i).^2);
    end
    
    disp(['error=',num2str(error), '  at time =',num2str(time)]);
    
    
    
    %if 1 && nargin==0 && mod(tstep,10)==0
    %   clf
    %    vv = Vp*p1;
    %    plot3(xp,yp,vv,'.');
    %    axis equal
    %    axis tight
    %    colorbar
    %    title(sprintf('time = %f',time))
    %    drawnow
    %end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
end



