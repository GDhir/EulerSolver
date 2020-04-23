function [ux,uy,w]=multisolver(X,Y,w,ux,uy,Nx,Ny,dx,dy,dt, cycle, XX, YY)


tol=10^-6;
omega=1.5;

%---------------Velocity Boundary conditions------------------%
for i=1:Nx
    for j=2:Ny-1
        for ii=2:Nx-1
            rx=XX(i,1)-XX(ii,j);
            ry=YY(i,1)-YY(ii,j);
            r=sqrt(rx^2+ry^2);
            ux(i,1)=ux(i,1)-ry*w(ii,j)*dx*dy/4/pi/(r^3);
            uy(i,1)=uy(i,1)+rx*w(ii,j)*dx*dy/4/pi/(r^3);
        end
    end
end

for j=1:Ny
    for jj=2:Ny-1
        for i=2:Nx-1
            rx=XX(1,j)-XX(i,jj);
            ry=YY(1,j)-YY(i,jj);
            r=sqrt(rx^2+ry^2);
            ux(1,j)=ux(1,j)-ry*w(i,jj)*dx*dy/4/pi/(r^3);
            uy(1,j)=uy(1,j)+rx*w(i,jj)*dx*dy/4/pi/(r^3);
        end
    end
end

for i=1:Nx
    for j=2:Ny-1
        for ii=2:Nx-1
            rx=XX(i,Ny)-XX(ii,j);
            ry=YY(i,Ny)-YY(ii,j);
            r=sqrt(rx^2+ry^2);
            ux(i,Ny)=ux(i,Ny)-ry*w(ii,j)*dx*dy/4/pi/(r^3);
            uy(i,Ny)=uy(i,Ny)+rx*w(ii,j)*dx*dy/4/pi/(r^3);
        end
    end
end

 for j=1:Ny
    for jj=2:Ny-1
        for i=2:Nx-1
            ry=XX(Nx,j)-XX(i,jj);
            rx=YY(Nx,j)-YY(i,jj);
            r=sqrt(rx^2+ry^2);
            ux(Nx,j)=ux(Nx,j)-ry*w(i,jj)*dx*dy/4/pi/(r^3);
            uy(Nx,j)=uy(Nx,j)+rx*w(i,jj)*dx*dy/4/pi/(r^3);
        end
    end
end

%----------------X Velocity field calculation-----------------%
error=10;
while error>tol
    uxo=ux;
    %SOR Relaxation
    for n=1:20
        uxo=ux;
        for j=2:Ny-1
            for i=2:Nx-1
                ux(i,j)=0.5*(dx^2*(w(i,j+1)-w(i,j-1))/2/dy+ux(i+1,j)+ux(i-1,j));
                ux(i,j)=omega*ux(i,j)+(1-omega)*uxo(i,j);
            end
        end
    end
    
    res=zeros(Nx,Ny);
    
    for j=2:Ny-1
        for i=2:Nx-1
            res(i,j)=-(w(i,j+1)-w(i,j-1))/2/dy-(ux(i+1,j)-2*ux(i,j)+ux(i-1,j))/dx^2;
        end
    end
    e(1)=mat2cell(ux,Nx,Ny);
    cycle=cycle+1;
    [res, e, Xp, Yp]=restrictionx(e, Nx, Ny, dx, res, X, Y, cycle);
    e=refinegrid(e,cycle,Xp,Yp);
    ux=cell2mat(e(1));
    error=sqrt(sum(sum((uxo-ux).^2))/Nx/Ny);
end


%----------------Y Velocity Field Calculation-----------------%
error=10;
while error>tol
    for n=1:20
        %SOR Relaxation
        uyo=uy;
        for j=2:Ny-1
            for i=2:Nx-1
                uy(i,j)=0.5*(-dy^2*(w(i+1,j)-w(i-1,j))/2/dx+uy(i,j+1)+uy(i,j-1));
                uy(i,j)=omega*uy(i,j)+(1-omega)*uyo(i,j);
            end
        end
    end
    res=zeros(Nx,Ny);
    for j=2:Ny-1
        for i=2:Nx-1
            res(i,j)=(w(i+1,j)-w(i-1,j))/2/dx-(uy(i,j+1)-2*uy(i,j)+uy(i,j-1))/dy^2;
        end
    end
     e(1)=mat2cell(uy,Nx,Ny);
    [res, e, Xp, Yp]=restrictionx(e, Nx, Ny, dy, res, X, Y, cycle);
    e=refinegrid(e,cycle,Xp,Yp);
    uy=cell2mat(e(1))
    
    error=sqrt(sum(sum((uyo-uy).^2))/Nx/Ny);
end    

%---------------Vorticity Field Calculation--------------------% 
wold=w;
for j=2:Ny-1
    for i=2:Nx-1
        w(i,j)=wold(i,j)-dt*(max(ux(i,j),0)*wold(i,j)+min(ux(i+1,j),0)*wold(i+1,j)-min(ux(i,j),0)*w(i,j)-max(ux(i-1,j),0)*w(i-1,j))/dx-...
            dt*(max(uy(i,j),0)*wold(i,j)+min(uy(i,j+1),0)*wold(i,j+1)-min(uy(i,j),0)*wold(i,j)-max(uy(i,j-1),0)*wold(i,j-1))/dy;
    end
end


