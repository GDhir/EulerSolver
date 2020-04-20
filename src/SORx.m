function u=SORx(Nx,Ny,u,h,res)

omega=1.9;
%SOR Relaxation
for n=1:20
    uo=u;
    for j=2:Ny-1
        for i=2:Nx-1
            u(i,j)=0.5*(u(i+1,j)+u(i-1,j)-h^2*res(i,j));
            u(i,j)=omega*u(i,j)+(1-omega)*uo(i,j);
        end
    end
end