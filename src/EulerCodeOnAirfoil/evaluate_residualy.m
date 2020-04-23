function res=evaluate_residualy(Nx,Ny,w,uy,dx,dy)

res=zeros(Nx,Ny);

    for j=2:Ny-1
        for i=2:Nx-1
            res(i,j)=(w(i+1,j)-w(i-1,j))/2/dx-(uy(i,j+1)-2*uy(i,j)+uy(i,j-1))/dy^2;
        end
    end