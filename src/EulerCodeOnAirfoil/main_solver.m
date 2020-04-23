function [q1, q2, q3, q4]=main_solver(vol, nx, ny, n_xi, n_eta, jtel , jteu, q1, q2, q3, q4, CFL, l, xc, yc)

global gamma rho_inf p_inf

Lp=zeros(4,4);
R=zeros(4,4);
Ln=zeros(4,4);
n=[1,0,-1,0];
m=[0,1,0,-1];
s=gamma-1;
error=10;
tol=10^-6;
dt=0.000000005;
dt1=10;
iter=0;

while error>tol
    q1o=q1;
    q2o=q2;
    q3o=q3;
    q4o=q4;
    for j=2:n_eta
        for i=2:n_xi
            Qc=[q1o(i,j), q2o(i,j), q3o(i,j), q4o(i,j)];
            rhoc=Qc(1);
            uc=Qc(2)/Qc(1);
            vc=Qc(3)/Qc(1);
            Ec=Qc(4);
            pc=(gamma-1)*Ec-rhoc*(uc^2+vc^2)*(gamma-1)/2;
            ac=sqrt(gamma*pc/rhoc);
            Hc=(Ec+pc)/rhoc;
            F=zeros(4,1);
            for k=1:4
                Q=[q1o(i+n(k),j+m(k)), q2o(i+n(k),j+m(k)), q3o(i+n(k),j+m(k)), q4o(i+n(k),j+m(k))];
                rho=Q(1);
                u=Q(2)/Q(1);
                v=Q(3)/Q(1);
                E=Q(4);
                p=(gamma-1)*E-rho*(u^2+v^2)*(gamma-1)/2;
                H=(E+p)/rho;
                vn=u*nx(k,i,j)+v*ny(k,i,j);
                a=sqrt(gamma*p/rho);
                Ln(1,1)=min(vn-a,0);
                Ln(2,2)=min(vn,0);
                Ln(3,3)=min(vn+a,0);
                Ln(4,4)=min(vn,0);
                R(:,1)=[1, u-a*nx(k,i,j), v-a*ny(k,i,j), H-vn*a];
                R(:,2)=[1, u, v, 0.5*(u^2+v^2)];
                R(:,3)=[1, u+a*nx(k,i,j), v+a*ny(k,i,j), H+vn*a];
                R(:,4)=[0,ny(k,i,j),-nx(k,i,j),u*ny(k,i,j)-v*nx(k,i,j)];
                L(:,1)=[(s*0.5*(u^2+v^2)+a*vn)/2/(a^2), (a^2-s*0.5*(u^2+v^2))/(a^2), (s*0.5*(u^2+v^2)-a*vn)/(2*a^2), v*nx(k,i,j)-u*ny(k,i,j)];
                L(:,2)=[(-s*u-a*nx(k,i,j))/(2*a^2), s*u/a^2, (-s*u+a*nx(k,i,j))/(2*a^2), ny(k,i,j)];
                L(:,3)=[(-s*v-a*ny(k,i,j))/(2*a^2), s*v/a^2, (-s*v+a*ny(k,i,j))/(2*a^2), -nx(k,i,j)];
                L(:,4)=[s/(2*a^2), -s/(a^2), s/(2*a^2), 0];
                Fnx=(R*(Ln*L))*Q';
                
                
                vn=uc*nx(k,i,j)+vc*ny(k,i,j);
                Lp(1,1)=max(vn-ac,0);
                Lp(2,2)=max(vn,0);
                Lp(3,3)=max(vn+ac,0);
                Lp(4,4)=max(vn,0);
                R(:,1)=[1, uc-ac*nx(k,i,j), vc-ac*ny(k,i,j), Hc-vn*ac];
                R(:,2)=[1, uc, vc, 0.5*(uc^2+vc^2)];
                R(:,3)=[1, uc+ac*nx(k,i,j), vc+ac*ny(k,i,j), Hc+vn*ac];
                R(:,4)=[0,ny(k,i,j),-nx(k,i,j),uc*ny(k,i,j)-vc*nx(k,i,j)];
                L(:,1)=[(s*0.5*(uc^2+vc^2)+ac*vn)/2/(ac^2), (ac^2-s*0.5*(uc^2+vc^2))/(ac^2), (s*0.5*(uc^2+vc^2)-ac*vn)/(2*ac^2), vc*nx(k,i,j)-uc*ny(k,i,j)];
                L(:,2)=[(-s*uc-ac*nx(k,i,j))/(2*ac^2), s*uc/ac^2, (-s*uc+ac*nx(k,i,j))/(2*ac^2), ny(k,i,j)];
                L(:,3)=[(-s*vc-ac*ny(k,i,j))/(2*ac^2), s*vc/ac^2, (-s*vc+ac*ny(k,i,j))/(2*ac^2), -nx(k,i,j)];
                L(:,4)=[s/(2*ac^2), -s/(ac^2), s/(2*ac^2), 0];
                Fpx=(R*(Lp*L))*Qc';
                F=F+(Fnx+Fpx)*l(k,i,j);
            end            
            q1(i,j)=q1o(i,j)-F(1)*dt/vol(i,j);
            q2(i,j)=q2o(i,j)-F(2)*dt/vol(i,j);
            q3(i,j)=q3o(i,j)-F(3)*dt/vol(i,j);
            q4(i,j)=q4o(i,j)-F(4)*dt/vol(i,j);
            Qc=[q1(i,j), q2(i,j), q3(i,j), q4(i,j)];
            rhoc=Qc(1);
            uc=Qc(2)/Qc(1);
            vc=Qc(3)/Qc(1);
            Ec=Qc(4);
            pc=(gamma-1)*Ec-rhoc*(uc^2+vc^2)*(gamma-1)/2;
            ac=sqrt(gamma*pc/rhoc);
            %vxi=uc*nx(3,i,j)+vc*ny(3,i,j);
            %veta=uc*nx(4,i,j)+vc*ny(4,i,j);
            dt1=min(dt1,CFL*min(l(:,i,j))/max(abs(uc)+ac,abs(vc)+ac));
        end
    end
    dt=dt1;
    %--------------Ghost Cells at solid boundaries--------------%
    for i=jtel+1:jteu
        q1(i,1)=q1(i,2);
        q2(i,1)=q2(i,2)*(ny(4,i,2)^2-nx(4,i,2)^2)-2*q3(i,2)*nx(4,i,2)*ny(4,i,2);
        q3(i,1)=q3(i,2)*(nx(4,i,2)^2-ny(4,i,2)^2)-2*q2(i,2)*nx(4,i,2)*ny(4,i,2);
        q4(i,1)=q4(i,2);
    end
    
    q1(2:jtel,1)=q1(n_xi:-1:jteu+1,2);
    q2(2:jtel,1)=q2(n_xi:-1:jteu+1,2);
    q3(2:jtel,1)=q3(n_xi:-1:jteu+1,2);
    q4(2:jtel,1)=q4(n_xi:-1:jteu+1,2);
    
    q1(jteu+1:n_xi,1)=q1(jtel:-1:2,2);
    q2(jteu+1:n_xi,1)=q2(jtel:-1:2,2);
    q3(jteu+1:n_xi,1)=q3(jtel:-1:2,2);
    q4(jteu+1:n_xi,1)=q4(jtel:-1:2,2);
    
    %{
    pin=(gamma-1)*q4(jtel+1:jteu, n_eta)-(q2(jtel+1:jteu, n_eta).^2+q3(jtel+1:jteu, n_eta).^2)*(gamma-1)/2./q1(jtel+1:jteu, n_eta);
    q4(jtel+1:jteu, n_eta+1)=pin/(gamma-1)+(q2(jtel+1:jteu, n_eta+1).^2+q3(jtel+1:jteu, n_eta+1).^2)*(gamma-1)/2./q1(jtel+1:jteu, n_eta+1);
    
    q1(jteu+1:n_xi, n_eta+1)=q1(jteu+1:n_xi, n_eta);
    q2(jteu+1:n_xi, n_eta+1)=q2(jteu+1:n_xi, n_eta);
    q3(jteu+1:n_xi, n_eta+1)=q3(jteu+1:n_xi, n_eta);
    q4(jteu+1:n_xi, n_eta+1)=p_inf/(gamma-1)+(q2(jteu+1:n_xi, n_eta+1).^2+q3(jteu+1:n_xi, n_eta+1).^2)*(gamma-1)/2./q1(jteu+1:n_xi, n_eta+1);
    %q4(jteu+1:n_xi, n_eta+1)=q4(jteu+1:n_xi, n_eta);
    
    q1(1:jtel, n_eta+1)=q1(1:jtel, n_eta);
    q2(1:jtel, n_eta+1)=q2(1:jtel, n_eta);
    q3(1:jtel, n_eta+1)=q3(1:jtel, n_eta);
    q4(1:jtel, n_eta+1)=p_inf/(gamma-1)+(q2(1:jtel, n_eta+1).^2+q3(1:jtel, n_eta+1).^2)*(gamma-1)/2./q1(1:jtel, n_eta+1);
    %q4(1:jtel, n_eta+1)=q4(1:jtel, n_eta);
    
    q1(1, 2:n_eta)=q1(2, 2:n_eta);
    q2(1, 2:n_eta)=q2(2, 2:n_eta);
    q3(1, 2:n_eta)=q3(2, 2:n_eta);
    q4(1, 2:n_eta)=p_inf/(gamma-1)+(q2(1, 2:n_eta).^2+q3(1, 2:n_eta).^2)*(gamma-1)/2./q1(1, 2:n_eta);
    %q4(1, 2:n_eta)=q4(2, 2:n_eta);
    
    q1(n_xi+1, 2:n_eta)=q1(n_xi, 2:n_eta);
    q2(n_xi+1, 2:n_eta)=q2(n_xi, 2:n_eta);
    q3(n_xi+1, 2:n_eta)=q3(n_xi, 2:n_eta);
    q4(n_xi+1, 2:n_eta)=p_inf/(gamma-1)+(q2(n_xi+1, 2:n_eta).^2+q3(n_xi+1, 2:n_eta).^2)*(gamma-1)/2./q1(n_xi+1, 2:n_eta);
    %q4(n_xi+1, 2:n_eta)=q4(n_xi, 2:n_eta);
    %}
    error=sqrt(sum(sum((q2-q2o).^2))/n_eta/n_xi)
    iter=iter+1;
    if(rem(iter,100)==0)
        plots_hw8(q1,q2,q3,q4, n_xi, n_eta, xc, yc, jteu, iter)
        error
        pause(25)
    end
end
         