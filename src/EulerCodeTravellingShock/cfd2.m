clc
clear

Nx1=10;
Ny1=20;
x0=0;
y0=-1;
y1=1;
global gamma
gamma=1.4;
CFL=0.01;

[x, y, Nx, Ny]=gridgen(Nx1, Ny1, x0, y0, y1);
[nx,ny,l,vol,xc,yc,P,facex,facey]=geometry(x,y, Nx, Ny);
Q=init(Nx,Ny,Nx1);
[Q,error]=mainsolver2(Q,nx,ny,xc,yc,Nx,Ny,vol,l,P,CFL,facex,facey);


figure(1)
plot(x,y,'.')
line(x,y)
line(x',y')
axis equal

figure(2)
contourf(xc(2:Nx,2:Ny),yc(2:Nx,2:Ny),reshape(Q(1,2:Nx,2:Ny),Nx-1,Ny-1))

figure(3)
loglog(error)