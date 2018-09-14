%clear
lambda = 1;
mu = 1;

syms k b1p b2p b1s b2s x y t B1 B2 B3 B4 w

 

u1a = (1i*k*B1*exp(-k*b1p*y)*exp(1i*(k*x-w*t)));

u2a = (-k*b1p*B1*exp(-k*b1p*y)*exp(1i*(k*x-w*t)));

u1e = (1i*k*B2*exp(k*b2p*y) - k*b2s*B3*exp(k*b2s*y)) * exp(1i*(k*x-w*t));

u2e = (k*b2p*B2*exp(k*b2p*y) + 1i*k*B3*exp(k*b2s*y)) * exp(1i*(k*x-w*t));

% take the real part of displacement

u1a_R = real(u1a);

u2a_R = real(u2a) ;

u1e_R = real(u1e);

u2e_R = real(u2e);

% velocity

v1a = simplify(diff(u1a_R,t));

v2a = simplify(diff(u2a_R,t));

v1b = simplify(diff(u1e_R,t));

v2b = simplify(diff(u2e_R,t));

 

% strain = .5*(grad(U)+grad(U)^T)

u1ax = simplify(diff(u1a_R,x));

u2ay = simplify(diff(u2a_R,y));

u12axy = 0.5*(simplify(diff(u1a_R,y) + diff(u2a_R,x)));

 
u1bx = simplify(diff(u1e_R,x));

u2by = simplify(diff(u2e_R,y));

u12bxy = 0.5*(simplify(diff(u1e_R,y) + diff(u2e_R,x)));

% stress = CE

s1ax = lambda*(u1ax + u2ay)

s2ay = lambda*(u1ax + u2ay)

s12axy = 0


s1bx = lambda*(u1bx + u2by) + 2*mu*u1bx

s2by = lambda*(u1bx + u2by) + 2*mu*u2by

s12bxy = 2*mu*u12bxy

% check source f
fex = simplify(simplify(diff(v1b,t))-simplify(diff(s1bx,x))-simplify(diff(s12bxy,y)))
fey = simplify(simplify(diff(v2b,t))-simplify(diff(s2by,y))-simplify(diff(s12bxy,x)))

fax = simplify(diff(v1a,t))-simplify(diff(s1ax,x))-simplify(diff(s12axy,y))
fay = simplify(diff(v2a,t))-simplify(diff(s2ay,y))-simplify(diff(s12axy,x))

