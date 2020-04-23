function Q=init(Nx,Ny,Nx1)


Q=zeros(4,Nx+1,Ny+1);

for j=1:Ny+1
    for i=1:Nx1
        up(2)=1.824;
        up(3)=0;
        up(4)=2.222;
        up(1)=3.783;        
        Q(:,i,j)=up2Q(up);              
    end
end

for j=1:Ny+1
    for i=Nx1+1:Nx+1
        up(2)=0;
        up(3)=0;
        up(4)=1;
        up(1)=1.4;        
        Q(:,i,j)=up2Q(up);            
    end
end

