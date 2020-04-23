function [res1, e, Xp, Yp]=restrictionx(e, Nx1, Ny1, dx, res1, X, Y, cycle)


Xp(1)=mat2cell(X,Nx1,Ny1);
Yp(1)=mat2cell(Y,Nx1,Ny1);
X1=X;
Y1=Y;
for level=2:cycle

    % Make the grid coarse
    Nx2=Nx1;
    Ny2=Ny1;
    dx1=(2^(level-1))*dx;
    Nx1=ceil(Nx1/2);
    Ny1=ceil(Ny1/2);
    resold=res1;
    res2=zeros(Nx1,Ny1);
    X2=zeros(Nx1,Ny1);
    Y2=zeros(Nx1,Ny1);
    if rem(Nx2,2)==0&&rem(Ny2,2)==0
        for j=1:Ny1
            for i=1:Nx1
                res2(i,j)=(resold(2*i-1,2*j-1)+resold(2*i-1,2*j)+resold(2*i,2*j-1)+resold(2*i,2*j))/4;
                X2(i,j)=(X1(2*i-1,2*j-1)+X1(2*i-1,2*j)+X1(2*i,2*j-1)+X1(2*i,2*j))/4;
                Y2(i,j)=(Y1(2*i-1,2*j-1)+Y1(2*i-1,2*j)+Y1(2*i,2*j-1)+Y1(2*i,2*j))/4;
            end
        end
    elseif rem(Nx2,2)==0
        for j=1:Ny1
            for i=1:Nx1
                res2(i,j)=(resold(2*i-1,2*j-1)+resold(2*i,2*j-1))/2;
                X2(i,j)=(X1(2*i-1,2*j-1)+X1(2*i,2*j-1))/2;
                Y2(i,j)=(Y1(2*i-1,2*j-1)+Y1(2*i,2*j-1))/2;
            end
        end     
    elseif rem(Ny2,2)==0
         for j=1:Ny1
            for i=1:Nx1
                res2(i,j)=(resold(2*i-1,2*j-1)+resold(2*i-1,2*j))/2;
                X2(i,j)=(X1(2*i-1,2*j-1)+X1(2*i-1,2*j))/2;
                Y2(i,j)=(Y1(2*i-1,2*j-1)+Y1(2*i-1,2*j))/2;
            end
         end
    else
        for j=1:Ny1
            for i=1:Nx1
                res2(i,j)=resold(2*i-1,2*j-1);
                X2(i,j)=X1(2*i-1,2*j-1);
                Y2(i,j)=Y1(2*i-1,2*j-1);
            end
         end
    end
    error=zeros(Nx1, Ny1);
    clear res1 X1 Y1
    X1=X2;
    Y1=Y2;
    % Relax the error on the coarser mesh
    error=SORx(Nx1,Ny1,error,dx1,res2);
    % Compute the new residual to be sent
    res1=evaluate_residualx(Nx1,Ny1,res2,error,dx1);
    e(level)=mat2cell(error,Nx1,Ny1);
    Xp(level)=mat2cell(X1,Nx1,Ny1);
    Yp(level)=mat2cell(Y1,Nx1,Ny1);
end 
