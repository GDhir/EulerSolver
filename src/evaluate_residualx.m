function res=evaluate_residualx(Nx,Ny,S,u,h)

res=zeros(Nx,Ny);

    for j=2:Ny-1
        for i=2:Nx-1
            res(i,j)=S(i,j)-(u(i+1,j)-2*u(i,j)+u(i-1,j))/h^2;
        end
    end