function Q=mainsolver(Q,nx,ny,xc,yc,Nx,Ny,vol,l,P,CFL,facex,facey)

global gamma
n=[0,1,0,-1];
m=[-1,0,1,0];
dt=0.005;
T=0.5;

t=0;
nn=0;
while t<T
    
    for i=2:Nx
        Q(1,i,1)=Q(1,i,2);
        Q(2,i,1)=Q(2,i,2);
        Q(3,i,1)=-Q(3,i,2);
        Q(4,i,1)=Q(4,i,2);
        
        Q(1,i,Ny+1)=Q(1,i,Ny);
        Q(2,i,Ny+1)=Q(2,i,Ny);
        Q(3,i,Ny+1)=-Q(3,i,Ny);
        Q(4,i,Ny+1)=Q(4,i,Ny);
    end
    
    for j=2:Ny
        Q(1,Nx+1,j)=Q(1,Nx,j);
        Q(2,Nx+1,j)=Q(2,Nx,j)*(ny(1,Nx-1,j-1)^2-nx(1,Nx-1,j-1)^2)-2*Q(3,Nx,j)*nx(1,Nx-1,j-1)*ny(1,Nx-1,j-1);
        Q(3,Nx+1,j)=Q(3,Nx,j)*(nx(1,Nx-1,j-1)^2-ny(1,Nx-1,j-1)^2)-2*Q(1,Nx-1,j)*nx(1,Nx-1,j-1)*ny(1,Nx-1,j-1);
        Q(4,Nx+1,j)=Q(4,Nx,j);
        
        Q(1,1,j)=3.783;
        Q(2,1,j)=Q(2,2,j);
        Q(3,1,j)=Q(3,2,j);
        Q(4,1,j)=2.222/(gamma-1)+3.783*(Q(2,2,j)^2+Q(3,2,j)^2)/(Q(1,2,j)^2)/2;
    end
    
    
        
    for j=2:Ny
        for i=2:Nx
            QL=[Q(1,i,j), Q(2,i,j), Q(3,i,j), Q(4,i,j)];
            flux=0;
            for k=1:4
                QR=[Q(1,i+n(k),j+m(k)), Q(2,i+n(k),j+m(k)), Q(3,i+n(k),j+m(k)), Q(4,i+n(k),j+m(k))];
                flux=flux+roeflux(QL, QR, nx(k,i-1,j-1), ny(k,i-1,j-1))*l(k,i-1,j-1);
            end
            Q(1,i,j)=Q(1,i,j)-flux(1)*dt/vol(i-1,j-1);
            Q(2,i,j)=Q(2,i,j)-flux(2)*dt/vol(i-1,j-1);
            Q(3,i,j)=Q(3,i,j)-flux(3)*dt/vol(i-1,j-1);
            Q(4,i,j)=Q(4,i,j)-flux(4)*dt/vol(i-1,j-1);
            up=Q2up(Q(:,i,j));
            a=sqrt(gamma*up(4)/up(1));
            l1=abs(up(2)-a);
            l2=abs(up(2));
            l3=abs(up(3));
            l4=abs(up(3)-a);
            s=max(l1,max(l2,max(l3,l4)));
            d=2*vol(i-1,j-1)/P(i-1,j-1);
            dt=min(dt,CFL*d/s);
            
        end
    end
   
    t=t+dt;
    nn=nn+1;
    if rem(nn,5)==0
        figure(2)
        contourf(xc(2:Nx,2:Ny),yc(2:Nx,2:Ny),reshape(Q(1,2:Nx,2:Ny),Nx-1,Ny-1))
        pause(2)
    end
end

