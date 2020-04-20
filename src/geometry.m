function [vol, nx, ny, xc, yc, l]=geometry(J, xix, xiy, etax, etay, n_xi, n_eta, x, y)

vol=zeros(n_xi+1,n_eta+1);
xc=zeros(n_xi+1,n_eta+1);
yc=zeros(n_xi+1,n_eta+1);
nx=zeros(4,n_xi+1,n_eta+1);
ny=zeros(4,n_xi+1,n_eta+1);
l=zeros(4,n_xi+1,n_eta+1);


for j=2:n_eta
    for i=2:n_xi
        vol(i,j)=(1/J(i-1,j-1)+1/J(i-1,j)+1/J(i,j-1)+1/J(i,j))/4;
        xc(i,j)=(x(i-1,j-1)+x(i-1,j)+x(i,j-1)+x(i,j))/4;
        yc(i,j)=(y(i-1,j-1)+y(i-1,j)+y(i,j-1)+y(i,j))/4;
        nx1(4)=-(etax(i-1,j-1)+etax(i,j-1))/2;
        nx1(2)=(etax(i-1,j)+etax(i,j))/2;
        nx1(3)=-(xix(i-1,j-1)+xix(i-1,j))/2;
        nx1(1)=(xix(i,j-1)+xix(i,j))/2;
        
        ny1(4)=-(etay(i-1,j-1)+etay(i,j-1))/2;
        ny1(2)=(etay(i-1,j)+etay(i,j))/2;
        ny1(3)=-(xiy(i-1,j-1)+xiy(i-1,j))/2;
        ny1(1)=(xiy(i,j-1)+xiy(i,j))/2;
        
        nx(4,i,j)=nx1(4)/sqrt(nx1(4)^2+ny1(4)^2);
        nx(2,i,j)=nx1(2)/sqrt(nx1(2)^2+ny1(2)^2);
        nx(3,i,j)=nx1(3)/sqrt(nx1(3)^2+ny1(3)^2);
        nx(1,i,j)=nx1(1)/sqrt(nx1(1)^2+ny1(1)^2);

        ny(4,i,j)=ny1(4)/sqrt(nx1(4)^2+ny1(4)^2);
        ny(2,i,j)=ny1(2)/sqrt(nx1(2)^2+ny1(2)^2);
        ny(3,i,j)=ny1(3)/sqrt(nx1(3)^2+ny1(3)^2);
        ny(1,i,j)=ny1(1)/sqrt(nx1(1)^2+ny1(1)^2);
        
        dx(4)=(x(i-1,j-1)-x(i,j-1));
        dx(2)=(x(i-1,j)-x(i,j));
        dx(3)=(x(i-1,j-1)-x(i-1,j));
        dx(1)=(x(i,j-1)-x(i,j));
        
        dy(4)=(y(i-1,j-1)-y(i,j-1));
        dy(2)=(y(i-1,j)-y(i,j));
        dy(3)=(y(i-1,j-1)-y(i-1,j));
        dy(1)=(y(i,j-1)-y(i,j));
        
        l(4,i,j)=sqrt(dx(4)^2+dy(4)^2);
        l(2,i,j)=sqrt(dx(2)^2+dy(2)^2);
        l(3,i,j)=sqrt(dx(3)^2+dy(3)^2);
        l(1,i,j)=sqrt(dx(1)^2+dy(1)^2);

    end
end

