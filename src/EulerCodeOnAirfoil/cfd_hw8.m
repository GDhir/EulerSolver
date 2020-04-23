clc
clear
readgrid

[J, xix, xiy, etax, etay] = metrics(n_xi, n_eta, jtel, jteu, x, y);

global gamma alpha M p_inf rho_inf
gamma=1.4;
alpha=2*pi/180;
M=1.2;
p_inf=101325;
rho_inf=1.225;

CFL=0.9;
[vol, nx, ny, xc, yc, l]=geometry(J, xix, xiy, etax, etay, n_xi, n_eta, x, y);
[q1, q2, q3, q4]=IC(n_xi, n_eta, jtel , jteu, nx, ny);

[q1, q2, q3, q4]=main_solver(vol, nx, ny, n_xi, n_eta, jtel , jteu, q1, q2, q3, q4, CFL, l, xc, yc);


