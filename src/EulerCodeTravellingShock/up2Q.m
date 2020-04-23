function Q=up2Q(up)

global gamma
Q(1)=up(1);
Q(2)=up(1)*up(2);
Q(3)=up(1)*up(3);
Q(4)=(up(4)/(gamma-1)+up(1)*(up(2)^2+up(3)^2)/2); 

