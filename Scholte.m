% set parameter

global Ka Ke

w = 2;
c1p = 1;
c2p = sqrt(3);
c2s = 1;

c = 0.711;
B1 = -1i*0.35945;
B2 = -1i*0.81946;
B3 = 1;

k = w/c;

b1p = sqrt(1-c^2/c1p^2);
b2p = sqrt(1-c^2/c2p^2);
b2s = sqrt(1-c^2/c2s^2);


%global s1ax s2ay s12axy s1ex s2ey s12exy
%global f1a f2a f1e f2e

u1a = @(x,y,t) real(B1*k*w*exp(-b1p*k*y).*exp(k*x*1i-t*w*1i));
u2a = @(x,y,t) -imag(B1*b1p*k*w*exp(-b1p*k*y).*exp(k*x.*1i-t*w*1i));

u1e = @(x,y,t) imag(k*w*exp(k*x*1i-t*w*1i).*(B2*exp(b2p*k*y)*1i-B3*b2s*exp(b2s*k*y)));
u2e = @(x,y,t) imag(k*w*exp(k*x.*1i-t*w*1i).*(B2*b2p*exp(b2p*k*y)+B3*exp(b2s*k*y)*1i));

s1ax = @(x,y,t) - real(B1*k^2*exp(-b1p*k*y).*exp(k*x*1i - t*w*1i)) + real(B1*b1p^2*k^2*exp(-b1p*k*y).*exp(k*x*1i - t*w*1i));
s2ay = @(x,y,t) - real(B1*k^2*exp(-b1p*k*y).*exp(k*x*1i - t*w*1i)) + real(B1*b1p^2*k^2*exp(-b1p*k*y).*exp(k*x*1i - t*w*1i));
s12axy = @(x,y,t) 0;

s1ex = @(x,y,t) real(k^2*exp(k*x*1i - t*w*1i).*(B2*b2p^2*exp(b2p*k*y) + B3*b2s*exp(b2s*k*y)*1i)) - 3*imag(k^2*exp(k*x*1i - t*w*1i).*(B2*exp(b2p*k*y)*1i - B3*b2s*exp(b2s*k*y)));
s2ey = @(x,y,t) 3*real(k^2*exp(k*x*1i - t*w*1i).*(B2*b2p^2*exp(b2p*k*y) + B3*b2s*exp(b2s*k*y)*1i)) - imag(k^2*exp(k*x*1i - t*w*1i).*(B2*exp(b2p*k*y)*1i - B3*b2s*exp(b2s*k*y)));
s12exy = @(x,y,t) - real(exp(k*x*1i - t*w*1i).*(B3*b2s^2*k^2*exp(b2s*k*y) - B2*b2p*k^2*exp(b2p*k*y)*1i)) - imag(k^2*exp(k*x*1i - t*w*1i).*(B2*b2p*exp(b2p*k*y) + B3*exp(b2s*k*y)*1i));


f1a = @(x,y,t) - imag(B1*k^3*exp(-b1p*k*y).*exp(k*x*1i - t*w*1i)) + imag(B1*k*w^2*exp(-b1p*k*y).*exp(k*x*1i - t*w*1i)) + imag(B1*b1p^2*k^3*exp(-b1p*k*y).*exp(k*x*1i - t*w*1i));

f2a = @(x,y,t) - real(B1*b1p*k^3*exp(-b1p*k*y).*exp(k*x*1i - t*w*1i)) + real(B1*b1p^3*k^3*exp(-b1p*k*y).*exp(k*x*1i - t*w*1i)) + real(B1*b1p*k*w^2*exp(-b1p*k*y).*exp(k*x*1i - t*w*1i));

f1e = @(x,y,t) - 3*imag(B2*k^3*exp(b2p*k*y).*exp(k*x*1i)*exp(-t*w*1i)) - real(B3*b2s*k^3*exp(b2s*k*y).*exp(k*x*1i)*exp(-t*w*1i)) + imag(B2*k*w^2*exp(b2p*k*y).*exp(k*x*1i)*exp(-t*w*1i)) + 3*imag(B2*b2p^2*k^3*exp(b2p*k*y).*exp(k*x*1i)*exp(-t*w*1i)) + real(B3*b2s^3*k^3*exp(b2s*k*y).*exp(k*x*1i)*exp(-t*w*1i)) + real(B3*b2s*k*w^2*exp(b2s*k*y).*exp(k*x*1i)*exp(-t*w*1i));

f2e = @(x,y,t) - imag(B3*k^3*exp(b2s*k*y).*exp(k*x*1i)*exp(-t*w*1i)) + 3*real(B2*b2p*k^3*exp(b2p*k*y).*exp(k*x*1i)*exp(-t*w*1i)) + imag(B3*k*w^2*exp(b2s*k*y).*exp(k*x*1i)*exp(-t*w*1i)) + imag(B3*b2s^2*k^3*exp(b2s*k*y).*exp(k*x*1i)*exp(-t*w*1i)) - 3*real(B2*b2p^3*k^3*exp(b2p*k*y).*exp(k*x*1i)*exp(-t*w*1i)) - real(B2*b2p*k*w^2*exp(b2p*k*y).*exp(k*x*1i)*exp(-t*w*1i));








if 0

time = 0;

% compute time step size
CN = (N+1)*(N+2)/3; % guessing...
dt = 2/(max(2*mu(:)+lambda(:))*CN*max(Fscale(:)));
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

% outer time step loop
tstep = 0;


for tstep = 1:Nsteps
        
    time = tstep*dt;
    
    U{1}(:,Ka) = u1a(x(:,Ka),y(:,Ka),time);
    U{1}(:,Ke) = u1e(x(:,Ke),y(:,Ke),time);

    U{2}(:,Ka) = u2a(x(:,Ka),y(:,Ka),time);
    U{2}(:,Ke) = u2e(x(:,Ke),y(:,Ke),time);

    U{3}(:,Ka) = s1ax(x(:,Ka),y(:,Ka),time);
    U{3}(:,Ke) = s1ex(x(:,Ke),y(:,Ke),time);

    U{4}(:,Ka) = s2ay(x(:,Ka),y(:,Ka),time);
    U{4}(:,Ke) = s2ey(x(:,Ke),y(:,Ke),time);

    U{5}(:,Ka) = zeros(Np,length(Ka));
    U{5}(:,Ke) = s12exy(x(:,Ke),y(:,Ke),time);
    
  
    if 1 && (mod(tstep,10)==0 || tstep==Nsteps)
        clf
        
        pe = (U{3} + U{4})/2; % trace(S)
        % pe = U{3}; % trace(S)
        pa = U{3}; % trace(S)
        p(:,Ke) = pe(:,Ke);
        p(:,Ka) = pa(:,Ka);
        vv = Vp*p;
        color_line3(xp,yp,vv,vv,'.');
        axis tight
        title(sprintf('time = %f',time));
        colorbar;
        
        drawnow        
        
    end
    
    if mod(tstep,100)==0
        disp(sprintf('On timestep %d out of %d\n',tstep,round(FinalTime/dt)))
    end
end
end