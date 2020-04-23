function up=Q2up(Q)

global gamma
up(1)=Q(1);
up(2)=Q(2)/Q(1);
up(3)=Q(3)/Q(1);
up(4)=(gamma-1)*Q(4)-up(1)*(up(2)^2+up(3)^2)*(gamma-1)/2;