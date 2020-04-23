function flux=rotatedroeflux(uL, uR, nx, ny)


global gamma
eps=10^-5;
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
pR = (gamma-1)*(uR(4) - 0.5*rhoR*(vxR*vxR+vyR*vyR));
aR = (gamma*pR/rhoR)^0.5;
HR = (uR(4) + pR)/rhoR;

fL(1) = rhoL*vnL;
fL(2) = rhoL*vnL * vxL + pL*nx;
fL(3) = rhoL*vnL * vyL + pL*ny;
fL(4) = rhoL*vnL *  HL;

fR(1) = rhoR*vnR;
fR(2) = rhoR*vnR * vxR + pR*nx;
fR(3) = rhoR*vnR * vyR + pR*ny;
fR(4) = rhoR*vnR *  HR;

%Define n1 and n2, and compute alpha1 and alpha2: (4.2) in the original paper.
%(NB: n1 and n2 may need to be frozen at some point during 
%     a steady calculation to fully make it converge. For time-accurate 
%     calculation, this is fine.)
% NB: For a boundary face, set (nx2,ny2)=(nx,ny), (nx1,ny1)=(-ny,nx).

    abs_dq = sqrt( (vxR-vxL)^2+(vyR-vyL)^2 );
  if  abs_dq > eps
       nx1 = (vxR-vxL)/abs_dq;
       ny1 = (vyR-vyL)/abs_dq;
  else
    nx1 = -ny ;
    ny1 =  nx;
  end
    alpha1 = nx * nx1 + ny * ny1 ;
%   To make alpha1 always positive.
      temp = sign(alpha1);
      if temp==0
          temp=1;
      end
       nx1 = temp * nx1;
       ny1 = temp * ny1;
    alpha1 = temp * alpha1;

% Take n2 as perpendicular to n1.
       nx2 = -ny1;
       ny2 =  nx1;
    alpha2 = nx * nx2 + ny * ny2;
%   To make alpha2 always positive.
      temp = sign(alpha2);
      if temp==0
          temp=1;
      end
       nx2 = temp * nx2;
       ny2 = temp * ny2;
    alpha2 = temp * alpha2;

%Now we are going to compute the Roe flux with n2 as the normal
%and n1 as the tagent vector, with modified wave speeds (5.12)

%Compute the Roe Averages
     RT = sqrt(rhoR/rhoL);
    rho = RT*rhoL;
     vx = (vxL+RT*vxR)/(1+RT);
     vy = (vyL+RT*vyR)/(1+RT);
      H = ( HL+RT* HR)/(1+RT);
      a = sqrt( (gamma-1)*(H-0.5*(vx*vx+vy*vy)) );
     vn = vx*nx2+vy*ny2;
     vt = vx*nx1+vy*ny1;

%Wave Strengths (remember that n2 is the normal and n1 is the tangent.)
    vnL = vxL*nx2 + vyL*ny2;
    vnR = vxR*nx2 + vyR*ny2;
    vtL = vxL*nx1 + vyL*ny1;
    vtR = vxR*nx1 + vyR*ny1;

   drho = rhoR - rhoL ;
     dp =   pR - pL;
    dvn =  vnR - vnL;
    dvt =  vtR - vtL;

  dV(1) = (dp - rho*a*dvn )/(2*a*a);
  dV(2) =  rho*dvt/a;
  dV(3) =  drho - dp/(a*a);
  dV(4) = (dp + rho*a*dvn )/(2*a*a);

%Wave Speeds for Roe flux part.
    ws(1) = vn-a;
    ws(2) = vn;
    ws(3) = vn;
    ws(4) = vn+a;
  abs_ws  = abs(ws);

%Harten's Entropy Fix JCP(1983), 49, pp357-393:
%only for the nonlinear fields.
  dws(1) = 1/5;
   if abs_ws(1)<dws(1) 
       abs_ws(1) = 0.5*(abs_ws(1)*abs_ws(1)/dws(1)+dws(1));
   end
       dws(4) = 1/5;
   if abs_ws(4)<dws(4) 
       abs_ws(4) = 0.5*(abs_ws(4)*abs_ws(4)/dws(4)+dws(4));
   end

%HLL wave speeds, evaluated with [nx1,ny1] (=tangent wrt n2).
   SRp = max( 0, max(vtR + aR, vt + a));
   SLm = min( 0, min(vtL - aL, vt - a));
%Modified wave speeds for the Rotated-RHLL flux: (5.12) in the original paper.
   ws = alpha2*abs_ws - ( alpha2*(SRp+SLm)*ws + 2*alpha1*SRp*SLm )/ (SRp-SLm);

%Right Eigenvectors: with n2 as normal and n1 as tangent.
  tx = nx1;
  ty = ny1;

  Rv(1,1) = 1   ; 
  Rv(2,1) = vx - a*nx2;
  Rv(3,1) = vy - a*ny2;
  Rv(4,1) =  H - vn*a;

  Rv(1,2) = 0;
  Rv(2,2) = a*tx;
  Rv(3,2) = a*ty;
  Rv(4,2) = a*vt;

  Rv(1,3) = 1;
  Rv(2,3) = vx;
  Rv(3,3) = vy ;
  Rv(4,3) = 0.5*(vx*vx+vy*vy);

  Rv(1,4) = 1;
  Rv(2,4) = vx + a*nx2;
  Rv(3,4) = vy + a*ny2;
  Rv(4,4) =  H + vn*a;

%Dissipation Term: Roe dissipation with the modified wave speeds.
  diss =[0,0,0,0];
  for i=1:4
      for j=1:4
        diss(i) = diss(i) + ws(j)*dV(j)*Rv(i,j);
      end 
  end 

%Compute the Rotated-RHLL flux.
  flux = (SRp*fL - SLm*fR)/(SRp-SLm) - 0.5*diss;

 
