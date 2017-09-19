function [rhsp, rhsu, rhsv] = acousticsRHS2D_WADG1(p,u,v)

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

rhsp = Pq*(Cq1.*(Vq*rhsp));
return;