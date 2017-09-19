function [rhsp, rhsu, rhsv] = acousticsRHS2D_adaptiveBBprojection(p,u,v, C_BB_proj,m)

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


pB = inv(T)*p; uB = inv(T)*u; vB = inv(T)*v;
%pr = Dr*p; ps = Ds*p;
pBr = 0.5*(-D0+D1)*pB; pBs = 0.5*(-D0+D2)*pB;

dpdx = rx.*(T*pBr) + sx.*(T*pBs);
%dpdx=-0.5*p0.*rx-0.5*p0.*sx+0.5*p1.*rx+0.5*p2.*sx;
dpdy = ry.*(T*pBr) + sy.*(T*pBs);
%dpdy=-0.5*p0.*ry-0.5*p0.*sy+0.5*p1.*ry+0.5*p2.*sy;
%divU = Dr*(u.*rx + v.*ry) + Ds*(u.*sx + v.*sy);
%divU = rx.*(Dr*u) + ry.*(Dr*v) + sx.*(Ds*u) + sy.*(Ds*v)
%divU = (D0*(inv(T)*u).*(-0.5*rx-0.5*sx)+inv(T)*v.*(-0.5*ry-0.5*sy))+...
%    D1*(inv(T)*u.*(0.5*rx)+inv(T)*v.*(0.5*ry))+D2*(inv(T)*u.*(0.5*sx)+inv(T)*v.*(0.5*sy));

divU = ((-0.5)*D0*uB).*rx + (0.5*D1*uB).*rx + ((-0.5)*D0*uB).*sx + (0.5*D2*uB).*sx...
      +((-0.5)*D0*vB).*ry + (0.5*D1*vB).*ry + ((-0.5)*D0*vB).*sy + (0.5*D2*vB).*sy;


% compute right hand sides of the PDE's
% I need to change the matrix LIFT


%rhsp =  -divU + LIFT1*T*inv(V)*(Fscale.*fluxp)/2.0;
%rhsu =  -dpdx + LIFT1*T*inv(V)*(Fscale.*fluxu)/2.0;
%rhsv =  -dpdy + LIFT1*T*inv(V)*(Fscale.*fluxv)/2.0;


rhsp =  -divU; % Notice: this coefficient is based on BB basis %+ LIFT1*inv(T)*(Fscale.*fluxp)/2.0;
rhsu =  -dpdx; %+ LIFT1*inv(T)*(Fscale.*fluxu)/2.0;
rhsv =  -dpdy; %+ LIFT1*inv(T)*(Fscale.*fluxv)/2.0;

[T0,T1,T2] = Vandermonde_BB_surface(r,s);

% Here rhsp is BB_coefficient
rhsp = rhsp + (LIFT1)* [inv(T0)*(Fscale(1:5,:).*fluxp(1:5,:));inv(T1)*(Fscale(6:10,:).*fluxp(6:10,:));inv(T2)*(Fscale(11:15,:).*fluxp(11:15,:))]/2.0;;


rhsu = rhsu + (T*LIFT1)* [inv(T0)*(Fscale(1:5,:).*fluxu(1:5,:));inv(T1)*(Fscale(6:10,:).*fluxu(6:10,:));inv(T2)*(Fscale(11:15,:).*fluxu(11:15,:))]/2.0;;


rhsv = rhsv + (T*LIFT1)* [inv(T0)*(Fscale(1:5,:).*fluxv(1:5,:));inv(T1)*(Fscale(6:10,:).*fluxv(6:10,:));inv(T2)*(Fscale(11:15,:).*fluxv(11:15,:))]/2.0;;


wu = zeros(15,512);

for i =0:4
P{i+1} = projection_BB_2D(4,i);
end



for j=1:512
    mp = (m(j)+1)*(m(j)+2)/2;
    
    wu(:,j) = P{m(j)+1} * mult_BB_2D(C_BB_proj(1:mp,j),rhsp(:,j),m(j),4);

   % wu(:,j) = P{m(j)+1}*wu(:,j);
end

%transform wu into nodal basis

rhsp = T * wu;  


%rhsp = Pq*(C.*(Vq*rhsp));



return;

