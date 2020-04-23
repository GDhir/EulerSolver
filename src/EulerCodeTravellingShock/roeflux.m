function flux=roeflux(uL, uR, nx, ny)

global gamma
tx=-ny;
ty=nx;

%Left State
rhoL=uL(1);
vxL = uL(2)/uL(1);
vyL = uL(3)/uL(1);
vnL = vxL*nx+vyL*ny;
vtL = vxL*tx+vyL*ty;
pL = (gamma-1)*(uL(4) - 0.5*rhoL*(vxL*vxL+vyL*vyL));
aL = (gamma*pL/rhoL)^0.5;
HL = (uL(4) + pL)/rhoL;

%Right State
rhoR = uR(1);
vxR = uR(2)/uR(1);
vyR = uR(3)/uR(1);
vnR = vxR*nx+vyR*ny;
vtR = vxR*tx+vyR*ty;
pR = (gamma-1)*(uR(3) - 0.5*rhoR*(vxR*vxR+vyR*vyR));
aR = (gamma*pR/rhoR)^0.5;
HR = (uR(4) + pR)/rhoR;

%Compute Roe's averages
RT = (rhoR/rhoL)^0.5;
rho = RT*rhoL;
vx = (vxL+RT*vxR)/(1+RT);
vy = (vyL+RT*vyR)/(1+RT);
H = ( HL+RT* HR)/(1+RT);
a = ((gamma-1)*(H-0.5*(vx*vx+vy*vy)))^0.5;
vn = vx*nx+vy*ny;
vt = vx*tx+vy*ty;

%Wave Strengths
drho = rhoR-rhoL;
dp = pR-pL;
dvn = vnR - vnL;
dvt = vtR - vtL;

dV(1) = (dp - rho*a*dvn)/(2*a*a);
dV(2) = rho*dvt/a;
dV(3) =  drho - dp/(a*a);
dV(4) = (dp + rho*a*dvn)/(2*a*a);

%Wave Speed
ws(1) = abs(vn-a);
ws(2) = abs(vn);
ws(3) = abs(vn);
ws(4) = abs(vn+a);

%Right eigenvectors in a Matrix
Rv(1, 1) = 1;
Rv(2, 1) = vx - a*nx;
Rv(3, 1) = vy - a*ny;
Rv(4, 1) = H - vn*a;

Rv(1, 2) = 0;
Rv(2, 2) = a*tx;
Rv(3, 2) = a*ty;
Rv(4, 2) = vt*a;

Rv(1, 3) = 1;
Rv(2, 3) = vx;
Rv(3, 3) = vy;
Rv(4, 3) = 0.5*(vx*vx+vy*vy);

Rv(1, 4) = 1;
Rv(2, 4) = vx + a*nx;
Rv(3, 4) = vy + a*ny;
Rv(4, 4) = H + vn*a;

diss=zeros(1,4);
%Dissipation Term
for i=1:4
    for j=1:4
        diss(1,i) = diss(1,i) + ws(j)*dV(j)*Rv(i, j);
    end
end
%Physical Fluxes
fluxL(1) = rhoL*vnL;
fluxL(2) = rhoL*vnL*vxL + pL*nx;
fluxL(3) = rhoL*vnL*vyL + pL*ny;
fluxL(4) = rhoL*vnL*HL;

fluxR(1) = rhoR*vnR;
fluxR(2) = rhoR*vnR*vxR + pR*nx;
fluxR(3) = rhoR*vnR*vyR + pR*ny;
fluxR(4) = rhoR*vnR*HR;

flux = 0.5*(fluxL + fluxR - diss);