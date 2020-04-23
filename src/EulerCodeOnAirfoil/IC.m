function [q1, q2, q3, q4]=IC(n_xi, n_eta, jtel , jteu, nx, ny)

global gamma alpha M p_inf rho_inf  

a_inf=sqrt(gamma*p_inf/rho_inf);
u=M*a_inf*cos(alpha);
v=M*a_inf*sin(alpha);

%------- Flow Domain and Ghost Cells for outer boundaries--------%
q1(1:n_xi+1,2:n_eta+1)=rho_inf*ones(n_xi+1,n_eta);
q2(1:n_xi+1,2:n_eta+1)=rho_inf*u*ones(n_xi+1,n_eta);
q3(1:n_xi+1,2:n_eta+1)=rho_inf*v*ones(n_xi+1,n_eta);
q4(1:n_xi+1,2:n_eta+1)=(p_inf/(gamma-1)+rho_inf*(u^2+v^2)/2)*ones(n_xi+1,n_eta);


%-------Ghost Cells Wake Cut----------%
q1(1:jtel,1)=rho_inf*ones(jtel,1);
q2(1:jtel,1)=rho_inf*u*ones(jtel,1);
q3(1:jtel,1)=rho_inf*v*ones(jtel,1);
q4(1:jtel,1)=(p_inf/(gamma-1)+rho_inf*(u^2+v^2)/2)*ones(jtel,1);

q1(jteu+1:n_xi+1,1)=rho_inf*ones(n_xi-jteu+1,1);
q2(jteu+1:n_xi+1,1)=rho_inf*u*ones(n_xi-jteu+1,1);
q3(jteu+1:n_xi+1,1)=rho_inf*v*ones(n_xi-jteu+1,1);
q4(jteu+1:n_xi+1,1)=(p_inf/(gamma-1)+rho_inf*(u^2+v^2)/2)*ones(n_xi-jteu+1,1);

%-------Ghost Cells Solid Boundaries------%
for i=jtel+1:jteu
    q1(i,1)=q1(i,2);
    q2(i,1)=q2(i,2)*(ny(4,i,2)^2-nx(4,i,2)^2)-2*q3(i,2)*nx(4,i,2)*ny(4,i,2);
    q3(i,1)=q3(i,2)*(nx(4,i,2)^2-ny(4,i,2)^2)-2*q2(i,2)*nx(4,i,2)*ny(4,i,2);
    q4(i,1)=q4(i,2);
end
