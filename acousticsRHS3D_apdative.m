function [rhsp, rhsu1, rhsu2, rhsu3] = acousticsRHS3D_apdative(p,u1,u2,u3)

Globals3D;

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du1 = zeros(Nfp*Nfaces,K); du1(:) = u1(vmapP)-u1(vmapM);
du2 = zeros(Nfp*Nfaces,K); du2(:) = u2(vmapP)-u2(vmapM);
du3 = zeros(Nfp*Nfaces,K); du3(:) = u3(vmapP)-u3(vmapM);

% evaluate upwind fluxes
ndotdU = nx.*du1 + ny.*du2 + nz.*du3;

% Impose reflective boundary conditions (p+ = -p-)
ndotdU(mapB) = 0;
dp(mapB) = -2*p(vmapB);

tau = 1;
fluxp =  tau*dp - ndotdU;
fluxu1 =  (tau*ndotdU - dp).*nx;
fluxu2 =  (tau*ndotdU - dp).*ny;
fluxu3 =  (tau*ndotdU - dp).*nz;

pr = Dr*p; ps = Ds*p; pt = Dt*p;
dpdx = rx.*pr + sx.*ps +tx.*pt;
dpdy = ry.*pr + sy.*ps +ty.*pt;
dpdz = rz.*pr + sz.*ps +tz.*pt;
divU = Dr*(u1.*rx + u2.*ry + u3.*rz) + Ds*(u1.*sx + u2.*sy + u3.*sz)+ Dt*(u1.*tx + u2.*ty + u3.*tz);

% compute right hand sides of the PDE's
rhsp =  -divU + LIFT*(Fscale.*fluxp)/2.0;
rhsu1 =  -dpdx + LIFT*(Fscale.*fluxu1)/2.0;
rhsu2 =  -dpdy + LIFT*(Fscale.*fluxu2)/2.0;
rhsu3 =  -dpdz + LIFT*(Fscale.*fluxu3)/2.0;

rhsp_BB = inv(TB)*rhsp;

for i=1:K
    
    k=mD(i);
    
    if (k==0) 
        v = inv(TB)*rhsp(:,i);
        c = v(1);
        rhsp(:,i)=rhsp(:,i)*c;
    end
    
    if (k==1 || k==2)
        kp = (k+1)*(k+2)*(k+3)/6;
        w = inv(TB)*rhsp(:,i);
        w1 = mult_BB_3D(w,c_BB(1:kp,i),N,k);
        Q = BB_projection3D(N,k);
        rhsp(:,i) = TB*Q*w1;
    end
    
    if(k==3)
        rhsp(:,i) = Pq*(Cq(:,i).*(Vq*rhsp(:,i)));
    end
    
end


return;