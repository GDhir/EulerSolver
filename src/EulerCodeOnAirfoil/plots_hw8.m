function []=plots_hw8(q1,q2,q3,q4, n_xi, n_eta, xc, yc, jteu, iter)
global gamma rho_inf p_inf

M_inf=1.2;

a_inf=sqrt(gamma*p_inf/rho_inf);
u_inf=M_inf*a_inf;


rhoc=q1;
uc=q2./q1;
vc=q3./q1;
Ec=q4;
pc=(gamma-1)*Ec-rhoc.*(uc.^2+vc.^2)*(gamma-1)/2;
Cp=(pc-p_inf)./(0.5*rho_inf*u_inf^2);


ac=sqrt(gamma*pc./rhoc);

M=sqrt(uc.^2+vc.^2)./ac;
figure(5)
plot(xc(100:jteu,2),Cp(100:jteu,2))
filename1 = ['xcp ' num2str(iter) '.jpg'];
 saveas(gcf,filename1)

figure(6)
contourf(xc(2:n_xi,2:n_eta),yc(2:n_xi,2:n_eta),q1(2:n_xi,2:n_eta),200)
filename2 = ['rho' num2str(iter) '.jpg'];
 saveas(gcf,filename2)


figure(7)
contourf(xc(2:n_xi,2:n_eta),yc(2:n_xi,2:n_eta),M(2:n_xi,2:n_eta),200)
filename1 = ['M' num2str(iter) '.jpg'];
 saveas(gcf,filename1)

figure(8)
contourf(xc(2:n_xi,2:n_eta),yc(2:n_xi,2:n_eta),Cp(2:n_xi,2:n_eta),200)
filename1 = ['cp' num2str(iter) '.jpg'];
 saveas(gcf,filename1)