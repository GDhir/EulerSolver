clc
clear
% Mesh Dimensions
n_xi     = 201;
n_eta    = 41;

% Trailing edge index on lower surface
jtel = 31; 

% Trailing edge index on upper surface
jteu = n_xi-jtel+1;

% Read grid
    gridxy=load('airfoil.dat');
    counter =0;
    for j=1:n_eta
	for i=1:n_xi
	   counter=counter+1;
	   x(i,j)=gridxy(counter,1);	
	   y(i,j)=gridxy(counter,2);	
	end
    end

% Plot grid
    figure(1)
    plot(x,y,'.')
    line(x,y)
    line(x',y')
    axis equal

% Zoom
    figure(2)
    plot(x,y,'.')
    line(x,y)
    line(x',y')
    xlim([-1 1])
    ylim([-1 1])

% Plot airfoil surface
    figure(3)
    plot(x(jtel:jteu,1),y(jtel:jteu,1),'o-')
    axis equal
