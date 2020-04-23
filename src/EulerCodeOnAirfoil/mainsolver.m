function [ux,uy,w]=mainsolver(X,Y,w,ux,uy,Nx,Ny,dx,dy,dt)


tol=10^-6;
omega=1.9;
ux(:,1)=zeros(Nx,1);
ux(:,Ny)=zeros(Nx,1);
uy(:,1)=zeros(Nx,1);
uy(:,Ny)=zeros(Nx,1);
ux(1,:)=zeros(1,Ny);
ux(Nx,:)=zeros(1,Ny);
uy(1,:)=zeros(1,Ny);
uy(Nx,:)=zeros(1,Ny);
%---------------Velocity Boundary conditions------------------%
for i=1:Nx
    for j=2:Ny-1
        for ii=2:Nx-1
            rx=X(i,1)-X(ii,j);
            ry=Y(i,1)-Y(ii,j);
            r=sqrt(rx^2+ry^2);
            ux(i,1)=ux(i,1)-ry*w(ii,j)*dx*dy/2/pi/(r);
            uy(i,1)=uy(i,1)+rx*w(ii,j)*dx*dy/2/pi/(r);
        end
    end
end

for j=1:Ny
    for jj=2:Ny-1
        for i=2:Nx-1
            rx=X(1,j)-X(i,jj);
            ry=Y(1,j)-Y(i,jj);
            r=sqrt(rx^2+ry^2);
            ux(1,j)=ux(1,j)-ry*w(i,jj)*dx*dy/2/pi/(r);
            uy(1,j)=uy(1,j)+rx*w(i,jj)*dx*dy/2/pi/(r);
        end
    end
end

for i=1:Nx
    for j=2:Ny-1
        for ii=2:Nx-1
            rx=X(i,Ny)-X(ii,j);
            ry=Y(i,Ny)-Y(ii,j);
            r=sqrt(rx^2+ry^2);
            ux(i,Ny)=ux(i,Ny)-ry*w(ii,j)*dx*dy/2/pi/(r);
            uy(i,Ny)=uy(i,Ny)+rx*w(ii,j)*dx*dy/2/pi/(r);
        end
    end
end

 for j=1:Ny
    for jj=2:Ny-1
        for i=2:Nx-1
            rx=X(Nx,j)-X(i,jj);
            ry=Y(Nx,j)-Y(i,jj);
            r=sqrt(rx^2+ry^2);
            ux(Nx,j)=ux(Nx,j)-ry*w(i,jj)*dx*dy/2/pi/(r);
            uy(Nx,j)=uy(Nx,j)+rx*w(i,jj)*dx*dy/2/pi/(r);
        end
    end
end

%----------------X Velocity field calculation-----------------%
error=10;
while error>tol
    uxo=ux;
    for j=2:Ny-1
        for i=2:Nx-1
            ux(i,j)=0.5*(dx^2*(w(i,j+1)-w(i,j-1))/2/dy+ux(i+1,j)+ux(i-1,j));
            ux(i,j)=omega*ux(i,j)+(1-omega)*uxo(i,j);
        end
    end
    error=sqrt(sum(sum((uxo-ux).^2))/Nx/Ny);
end

%----------------Y Velocity Field Calculation-----------------%
error=10;
while error>tol
    uyo=uy;
    for j=2:Ny-1
        for i=2:Nx-1
            uy(i,j)=0.5*(-dy^2*(w(i+1,j)-w(i-1,j))/2/dx+uy(i,j+1)+uy(i,j-1));
            uy(i,j)=omega*uy(i,j)+(1-omega)*uyo(i,j);
        end
    end
    error=sqrt(sum(sum((uyo-uy).^2))/Nx/Ny);
end    

%---------------Vorticity Field Calculation--------------------% 
wold=w;
for j=2:Ny-1
    for i=2:Nx-1
        w(i,j)=wold(i,j)-dt*(max(ux(i,j),0)*wold(i,j)+min(ux(i+1,j),0)*wold(i+1,j)-min(ux(i,j),0)*wold(i,j)-max(ux(i-1,j),0)*wold(i-1,j))/dx-...
            dt*(max(uy(i,j),0)*wold(i,j)+min(uy(i,j+1),0)*wold(i,j+1)-min(uy(i,j),0)*wold(i,j)-max(uy(i,j-1),0)*wold(i,j-1))/dy;
    end
end


