function[J, xix, xiy, etax, etay] = metrics(n_xi, n_eta, jtel, jteu, x, y)

J=zeros(n_xi,n_eta);

xxi=zeros(n_xi,n_eta);
xeta=zeros(n_xi,n_eta);
yxi=zeros(n_xi,n_eta);
yeta=zeros(n_xi,n_eta);

% Inner points
 for j=2:n_eta-1
     for i=2:n_xi-1
        xxi(i,j)=(x(i+1,j)-x(i-1,j))/2;
        yxi(i,j)=(y(i+1,j)-y(i-1,j))/2;
        xeta(i,j)=(x(i,j+1)-x(i,j-1))/2;
        yeta(i,j)=(y(i,j+1)-y(i,j-1))/2;
        J(i,j)=1/(xxi(i,j)*yeta(i,j)-yxi(i,j)*xeta(i,j));
     end
 end
 clear i j;
 
 % Inner Wake Points 
 for i=2:jtel-1
        xxi(i,1)=(x(i+1,1)-x(i-1,1))/2;
        yxi(i,1)=(y(i+1,1)-y(i-1,1))/2;
        xeta(i,1)=(x(i,2)-x(jteu+jtel-i,2))/2;
        yeta(i,1)=(y(i,2)-y(jteu+jtel-i,2))/2; 
        J(i,1)=1/(xxi(i,1)*yeta(i,1)-yxi(i,1)*xeta(i,1));
 end
 
 clear i j;
 
 
 % Airfoil Surface
 for i=jtel:jteu
        xxi(i,1)=(x(i+1,1)-x(i-1,1))/2;
        yxi(i,1)=(y(i+1,1)-y(i-1,1))/2;
        xeta(i,1)=(-3*x(i,1)+4*x(i,2)-x(i,3))/2;
        yeta(i,1)=(-3*y(i,1)+4*y(i,2)-y(i,3))/2; 
        J(i,1)=1/(xxi(i,1)*yeta(i,1)-yxi(i,1)*xeta(i,1));
 end
 
 clear i j;
 
 % Inner Wake Points 
 for i=jteu+1:n_xi-1
        xxi(i,1)=(x(i+1,1)-x(i-1,1))/2;
        yxi(i,1)=(y(i+1,1)-y(i-1,1))/2;
        xeta(i,1)=(x(i,2)-x(jteu+jtel-i,2))/2;
        yeta(i,1)=(y(i,2)-y(jteu+jtel-i,2))/2; 
        J(i,1)=1/(xxi(i,1)*yeta(i,1)-yxi(i,1)*xeta(i,1));
 end
 
 clear i j;
 
 
 
 % Upper Boundary
 for i=2:n_xi-1
     xxi(i,n_eta)=(x(i+1,n_eta)-x(i-1,n_eta))/2;
     yxi(i,n_eta)=(y(i+1,n_eta)-y(i-1,n_eta))/2;
     xeta(i,n_eta)=(3*x(i,n_eta)-4*x(i,n_eta-1)+x(i,n_eta-2))/2;
     yeta(i,n_eta)=(3*y(i,n_eta)-4*y(i,n_eta-1)+y(i,n_eta-2))/2; 
     J(i,n_eta)=1/(xxi(i,n_eta)*yeta(i,n_eta)-yxi(i,n_eta)*xeta(i,n_eta));
 end
 
    clear i j;
 
 % Lower Side Boundary
 for j=2:n_eta-1
     xxi(1,j)=(-3*x(1,j)+4*x(2,j)-x(3,j))/2;
     yxi(1,j)=(-3*y(1,j)+4*y(2,j)-y(3,j))/2;
     xeta(1,j)=(x(1,j+1)-x(1,j-1))/2;
     yeta(1,j)=(y(1,j+1)-y(1,j-1))/2;
     J(1,j)=1/(xxi(1,j)*yeta(1,j)-yxi(1,j)*xeta(1,j));
 end
 
   clear i j;
 
 % Upper Side Boundary
 for j=2:n_eta-1
     xxi(n_xi,j)=(3*x(n_xi,j)-4*x(n_xi-1,j)+x(n_xi-2,j))/2;
     yxi(n_xi,j)=(3*y(n_xi,j)-4*y(n_xi-1,j)+y(n_xi-2,j))/2;
     xeta(n_xi,j)=(x(n_xi,j+1)-x(n_xi,j-1))/2;
     yeta(n_xi,j)=(y(n_xi,j+1)-y(n_xi,j-1))/2;
     J(n_xi,j)=1/(xxi(n_xi,j)*yeta(n_xi,j)-yxi(n_xi,j)*xeta(n_xi,j));
 end
   
    clear i j;
 
 % Corner Points
     xxi(1,n_eta)=(-3*x(1,n_eta)+4*x(2,n_eta)-x(3,n_eta))/2;
     yxi(1,n_eta)=(-3*y(1,n_eta)+4*y(2,n_eta)-y(3,n_eta))/2;
     xeta(1,n_eta)=(3*x(1,n_eta)-4*x(1,n_eta-1)+x(1,n_eta-2))/2;
     yeta(1,n_eta)=(3*y(1,n_eta)-4*y(1,n_eta-1)+y(1,n_eta-2))/2;
     J(1,n_eta)=1/(xxi(1,n_eta)*yeta(1,n_eta)-yxi(1,n_eta)*xeta(1,n_eta));
     
     xxi(n_xi,n_eta)=(3*x(n_xi,n_eta)-4*x(n_xi-1,n_eta)+x(n_xi-2,n_eta))/2;
     yxi(n_xi,n_eta)=(3*y(n_xi,n_eta)-4*y(n_xi-1,n_eta)+y(n_xi-2,n_eta))/2;
     xeta(n_xi,n_eta)=(3*x(n_xi,n_eta)-4*x(n_xi,n_eta-1)+x(n_xi,n_eta-2))/2;
     yeta(n_xi,n_eta)=(3*y(n_xi,n_eta)-4*y(n_xi,n_eta-1)+y(n_xi,n_eta-2))/2;
     J(n_xi,n_eta)=1/(xxi(n_xi,n_eta)*yeta(n_xi,n_eta)-yxi(n_xi,n_eta)*xeta(n_xi,n_eta));
     
  % Boundary Wake Points
     xxi(1,1)=(-3*x(1,1)+4*x(2,1)-x(3,1))/2;
     yxi(1,1)=(-3*y(1,1)+4*y(2,1)-y(3,1))/2;
     xeta(1,1)=(x(1,2)-x(jteu+jtel-1,2))/2;
     yeta(1,1)=(y(1,2)-y(jteu+jtel-1,2))/2; 
     J(1,1)=1/(xxi(1,1)*yeta(1,1)-yxi(1,1)*xeta(1,1));
     
     xxi(n_xi,1)=(3*x(n_xi,1)-4*x(n_xi-1,1)+x(n_xi-2,1))/2;
     yxi(n_xi,1)=(3*y(n_xi,1)-4*y(n_xi-1,1)+y(n_xi-2,1))/2;
     xeta(n_xi,1)=(x(n_xi,2)-x(jteu+jtel-n_xi,2))/2;
     yeta(n_xi,1)=(y(n_xi,2)-y(jteu+jtel-n_xi,2))/2; 
     J(n_xi,1)=1/(xxi(n_xi,1)*yeta(n_xi,1)-yxi(n_xi,1)*xeta(n_xi,1));
     
     
  xix=J.*yeta;
  etax=-J.*yxi;
  xiy=-J.*xeta;
  etay=J.*xxi;
  
  %contourf(x,y,etay)
        
 
 