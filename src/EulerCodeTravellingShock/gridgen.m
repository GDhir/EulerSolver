function [x, y, Nx, Ny]=gridgen(Nx1, Ny, x0, y0, y1)

Nx2=Nx1;
Nx=Nx1+Nx2;
y=zeros(Nx,Ny);
x=zeros(Nx,Ny);


dy=(y1-y0)/(Ny-1);


for j=1:Ny
    for i=1:Nx
        y(i,j)=y0+(j-1)*dy;
    end
end

x1dash=1;
dx1=(x1dash-x0)/(Nx1-1);
x(Nx,:)=2-(y(Nx,:).^2)/4;

dx0=dx1;
d2x(1,1:Ny)=2*((Nx2-1)*dx0-(x(Nx,:)-x1dash))/(Nx2-1)/(Nx2);

for j=1:Ny
    for i=1:Nx1
        x(i,j)=x0+(i-1)*dx1;
    end
end

for j=1:Ny
    for i=Nx1+1:Nx
        dx2=dx0-(i-Nx1-1)*d2x(1,j);
        x(i,j)=x(i-1,j)+dx2;
    end
end
