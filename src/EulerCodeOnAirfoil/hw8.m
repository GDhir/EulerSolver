clear
clc
Nx=100;
Ny=100;
xl=-20;
xr=20;
yl=-20;
yr=20;
x=linspace(xl,xr,Nx);
y=linspace(yl,yr,Ny);
dx=(xr-xl)/(Nx-1);
dy=(yr-yl)/(Ny-1);
[X,Y]=meshgrid(x,y);
X=X';
Y=Y';
w=zeros(Nx,Ny);
ux=zeros(Nx,Ny);
uy=zeros(Nx,Ny);
circ1=1;
circ2=0.1;
rc1=1;
rc2=0.5;
T=100;
CFL=0.9;


for j=1:Ny
    for i=1:Nx
        r1=sqrt((X(i,j)-(-5))^2+(Y(i,j)-0)^2);
        r2=sqrt((X(i,j)-(-1))^2+(Y(i,j)-0)^2);
        r3=sqrt((X(i,j)-(5))^2+(Y(i,j)-0)^2);
        r4=sqrt((X(i,j)-(1))^2+(Y(i,j)-0)^2);
        if r1<rc1
            w(i,j)=-circ1*exp(-(r1^2)/(rc1^2))/pi/(rc1^2);
        end
        
        if r2<rc2
            w(i,j)=-circ2*exp(-(r2^2)/(rc2^2))/pi/(rc2^2);
        end
        
        if r3<rc1
            w(i,j)=circ1*exp(-(r3^2)/(rc1^2))/pi/(rc1^2);
        end
        
        if r4<rc2
            w(i,j)=circ2*exp(-(r4^2)/(rc2^2))/pi/(rc2^2);
        end
    end
end
cycle=2;
t=0;
dt=10^-5;
while t<T 
    [ux,uy,w]=mainsolver(X,Y,w,ux,uy,Nx,Ny,dx,dy,dt);
    dt=CFL*min(dx/max(max(ux)),dy/max(max(uy)));
    t=t+dt
    V=ux.^2+uy.^2;
    contourf(X,Y,w,25)
    pause(5)
end