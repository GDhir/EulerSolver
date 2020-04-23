function [nx,ny,l,vol,xc, yc,P,facex,facey]=geometry(x,y, Nx, Ny)


nx=zeros(4,Nx-1,Ny-1);
ny=zeros(4,Nx-1,Ny-1);
l=zeros(4,Nx-1,Ny-1);
vol=zeros(Nx-1,Ny-1);
xc=zeros(Nx+1,Ny+1);
yc=zeros(Nx+1,Ny+1);
P=zeros(Nx-1,Ny-1);
facex=zeros(4,Nx+1,Ny+1);
facey=zeros(4,Nx+1,Ny+1);


for j=1:Ny-1
    for i=1:Nx-1
        dx(1)=x(i+1,j)-x(i,j);
        dx(2)=x(i+1,j+1)-x(i+1,j);
        dx(3)=x(i,j+1)-x(i+1,j+1);
        dx(4)=x(i,j)-x(i,j+1);
        xc(i+1,j+1)=(x(i,j)+x(i,j+1)+x(i+1,j)+x(i+1,j+1))/4;
        yc(i+1,j+1)=(y(i,j)+y(i,j+1)+y(i+1,j)+y(i+1,j+1))/4;
        
        facex(1,i+1,j+1)=(x(i,j)+x(i+1,j))/2;
        facex(2,i+1,j+1)=(x(i+1,j+1)+x(i+1,j))/2;
        facex(3,i+1,j+1)=(x(i,j+1)+x(i+1,j+1))/2;
        facex(4,i+1,j+1)=(x(i,j)+x(i,j+1))/2;
        
        facey(1,i+1,j+1)=(y(i,j)+y(i+1,j))/2;
        facey(2,i+1,j+1)=(y(i+1,j+1)+y(i+1,j))/2;
        facey(3,i+1,j+1)=(y(i,j+1)+y(i+1,j+1))/2;
        facey(4,i+1,j+1)=(y(i,j)+y(i,j+1))/2;
        
        
        dy(1)=y(i+1,j)-y(i,j);
        dy(2)=y(i+1,j+1)-y(i+1,j);
        dy(3)=y(i,j+1)-y(i+1,j+1);
        dy(4)=y(i,j)-y(i,j+1);
        P(i,j)=0;
        
        for k=1:4
            nx1(k)=dy(k);
            ny1(k)=-dx(k);
            ny(k,i,j)=ny1(k)/sqrt(nx1(k)^2+ny1(k)^2);
            nx(k,i,j)=nx1(k)/sqrt(nx1(k)^2+ny1(k)^2);
            l(k,i,j)=sqrt(dx(k)^2+dy(k)^2);
            P(i,j)=P(i,j)+l(k,i,j);
        end
        
        X=[x(i,j), x(i+1,j), x(i+1,j+1), x(i,j+1)];
        Y=[y(i,j), y(i+1,j), y(i+1,j+1), y(i,j+1)];
        
        vol(i,j)=polyarea(X,Y);
        if i==1
            xc(i,j+1)=xc(i+1,j+1)-dx(1);
            yc(i,j+1)=yc(i+1,j+1);
        elseif i==Nx-1
            xc(i+2,j+1)=xc(i+1,j+1)+dx(1);
            yc(i+2,j+1)=yc(i+1,j+1);
        elseif j==1
            xc(i+1,j)=xc(i+1,j+1);
            yc(i+1,j)=yc(i+1,j+1)-dy(2);
        elseif j==Ny-1
            xc(i+1,j+2)=xc(i+1,j+1);
            yc(i+1,j+2)=yc(i+1,j+1)+dy(2);
        end
        
    end
end
       

            