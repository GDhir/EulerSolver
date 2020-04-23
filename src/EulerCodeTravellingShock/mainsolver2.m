function [Q,error]=mainsolver2(Q,nx,ny,xc,yc,Nx,Ny,vol,l,P,CFL,facex,facey)

global gamma
n=[0,1,0,-1];
m=[-1,0,1,0];
ll=[3,4,1,2];
dt=0.0005;
T=0.05;
epsdot=0.05;
t=0;
nn=0;

% LSQ
a1=zeros(Nx+1,Ny+1);
b=zeros(Nx+1,Ny+1);
c=zeros(Nx+1,Ny+1);
dxik=zeros(4,Nx+1,Ny+1);
dyik=zeros(4,Nx+1,Ny+1);
wik=zeros(4,Nx+1,Ny+1);
phi=zeros(4,Nx+1,Ny+1);
for j=2:Ny
    for i=2:Nx
        for k=1:4
            dxik(k,i,j)=xc(i+n(k),j+m(k))-xc(i,j);
            dyik(k,i,j)=yc(i+n(k),j+m(k))-yc(i,j);
            wik(k,i,j)=1\sqrt(dxik(k,i,j)^2+dyik(k,i,j)^2);
            a1(i,j)=a1(i,j)+(wik(k,i,j)^2)*dxik(k,i,j)^2;
            c(i,j)=c(i,j)+(wik(k,i,j)^2)*dyik(k,i,j)^2;
            b(i,j)=b(i,j)+(wik(k,i,j)^2)*dxik(k,i,j)*dyik(k,i,j);
        end
    end
end



while t<T
    Qold=Q;
    upx=zeros(4,Nx+1,Ny+1);
    upy=zeros(4,Nx+1,Ny+1);
     
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
    
    % Find Global max and min values
    umaxg=[0,0,0,0];
    uming=[0,0,0,0];
    for j=2:Ny
        for i=2:Nx
            Qc=[Q(1,i,j), Q(2,i,j), Q(3,i,j), Q(4,i,j)];            
            uc=Q2up(Qc);
            umaxg(1)=max(umaxg(1),uc(1));
            umaxg(2)=max(umaxg(2),uc(2));
            umaxg(3)=max(umaxg(3),uc(3));
            umaxg(4)=max(umaxg(4),uc(4));
            
            uming(1)=min(uming(1),uc(1));
            uming(2)=min(uming(2),uc(2));
            uming(3)=min(uming(3),uc(3));
            uming(4)=min(uming(4),uc(4));
        end
    end
    
    % Find Gradients and local max and min values
    umax=zeros(4,Nx+1,Ny+1);
    umin=zeros(4,Nx+1,Ny+1);
    for j=2:Ny
        for i=2:Nx
            Qc=[Q(1,i,j), Q(2,i,j), Q(3,i,j), Q(4,i,j)];
            
            uc=Q2up(Qc);
            umax(:,i,j)=uc;
            umin(:,i,j)=uc;
            di=[0,0,0,0];
            ei=[0,0,0,0];
            

            for k=1:4
                Qik=[Q(1,i+n(k),j+m(k)),Q(2,i+n(k),j+m(k)),Q(3,i+n(k),j+m(k)),Q(4,i+n(k),j+m(k))];
                
                uik=Q2up(Qik);
                umax(1,i,j)=max(uik(1),umax(1,i,j));
                umax(2,i,j)=max(uik(2),umax(2,i,j));
                umax(3,i,j)=max(uik(3),umax(3,i,j));
                umax(4,i,j)=max(uik(4),umax(4,i,j));
                
                umin(1,i,j)=min(uik(1),umin(1,i,j));
                umin(2,i,j)=min(uik(2),umin(2,i,j));
                umin(3,i,j)=min(uik(3),umin(3,i,j));
                umin(4,i,j)=min(uik(4),umin(4,i,j));
                
                duik=uik-uc;
                di=di+(wik(k,i,j)^2)*duik*dxik(k,i,j);
                ei=ei+(wik(k,i,j)^2)*duik*dyik(k,i,j);
                %QR=[Q(1,i+n(k),j+m(k)), Q(2,i+n(k),j+m(k)), Q(3,i+n(k),j+m(k)), Q(4,i+n(k),j+m(k))];
                
            end
            
            
            upx(:,i,j)=(c(i,j)*di-b(i,j)*ei)/(a1(i,j)*c(i,j)-b(i,j)^2);
            upy(1:4,i,j)=(-b(i,j)*di(1:4)+a1(i,j)*ei(1:4))/(a1(i,j)*c(i,j)-b(i,j)^2);
        end
    end
    
    %Limiting
    upll=zeros(4,4,Nx+1,Ny+1);
    for j=2:Ny
        for i=2:Nx
            Qc=[Q(1,i,j), Q(2,i,j), Q(3,i,j), Q(4,i,j)];
            ucl=Q2up(Qc); 
            phiold=phi(:,i,j);
            for k=1:4
                upl(1)=ucl(1)+upx(1,i,j)*(facex(k,i,j)-xc(i,j))+upy(1,i,j)*(facey(k,i,j)-yc(i,j));
                upl(2)=ucl(2)+upx(2,i,j)*(facex(k,i,j)-xc(i,j))+upy(2,i,j)*(facey(k,i,j)-yc(i,j));
                upl(3)=ucl(3)+upx(3,i,j)*(facex(k,i,j)-xc(i,j))+upy(3,i,j)*(facey(k,i,j)-yc(i,j));
                upl(4)=ucl(4)+upx(4,i,j)*(facex(k,i,j)-xc(i,j))+upy(4,i,j)*(facey(k,i,j)-yc(i,j));
                
                %QL=up2Q(upl);
                
                deltan=upl-ucl;
                deltap=[0,0,0,0];
                for t=1:4
                    if deltan(t)>0
                        deltap(t)=umax(t,i,j)-ucl(t);
                    else
                        deltap(t)=umin(t,i,j)-ucl(t);
                    end
                end                
                eps=epsdot*(umaxg-uming);
                
                for t=1:4
                    if deltan(t)==0
                        phi(t,i,j)=1;
                    else
                        phi(t,i,j)=((deltan(t)^2+eps(t)^2)*deltan(t)+2*(deltan(t)^2)*deltap(t))/(deltap(t)^2+2*deltan(t)^2+deltap(t)*deltan(t)+eps(t)^2)/deltan(t);
                    end
                    phi(t,i,j)=min(phiold(t),phi(t,i,j));
                end
                
            end
             
        end
    end
    
    
    % Find reconstructed values
     for j=2:Ny
        for i=2:Nx
             Qc=[Q(1,i,j), Q(2,i,j), Q(3,i,j), Q(4,i,j)];
             ucl=Q2up(Qc); 
             for k=1:4
                    upll(1,k,i,j)=ucl(1)+phi(1,i,j)*(upx(1,i,j)*(facex(k,i,j)-xc(i,j))+upy(1,i,j)*(facey(k,i,j)-yc(i,j)));
                    upll(2,k,i,j)=ucl(2)+phi(2,i,j)*(upx(2,i,j)*(facex(k,i,j)-xc(i,j))+upy(2,i,j)*(facey(k,i,j)-yc(i,j)));
                    upll(3,k,i,j)=ucl(3)+phi(3,i,j)*(upx(3,i,j)*(facex(k,i,j)-xc(i,j))+upy(3,i,j)*(facey(k,i,j)-yc(i,j)));
                    upll(4,k,i,j)=ucl(4)+phi(4,i,j)*(upx(4,i,j)*(facex(k,i,j)-xc(i,j))+upy(4,i,j)*(facey(k,i,j)-yc(i,j)));                
             end
        end
     end
     
     
     
     
     for j=2:Ny
        for i=2:Nx
            flux=0;
            for k=1:4
                uplll=upll(:,k,i,j);
                uprrr=upll(:,ll(k),i+n(k),j+m(k));
                QR=up2Q(uprrr);
                QL=up2Q(uplll);
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
    error(nn)=sqrt(sum(sum((Q(1,:,:)-Qold(1,:,:)).^2))/Nx/Ny);
   %{
    if rem(nn,5)==0
        figure(2)
        contourf(xc(2:Nx,2:Ny),yc(2:Nx,2:Ny),reshape(Q(1,2:Nx,2:Ny),Nx-1,Ny-1))
        pause(2)
    %}
end
