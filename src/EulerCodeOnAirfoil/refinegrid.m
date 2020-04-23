function e=refinegrid(e,cycle,Xp,Yp)


    for level=cycle-1:-1:1
        X1=cell2mat(Xp(level+1));
        Y1=cell2mat(Yp(level+1));
        X2=cell2mat(Xp(level));
        Y2=cell2mat(Yp(level));
        temp=cell2mat(e(level+1));
        temp=(interp2(X1,Y1,temp',X2,Y2))';
        temp1=cell2mat(e(level));
        temp2=temp1+temp; 
        e(level)=mat2cell(temp2,size(temp2,1),size(temp2,2));
        clear temp;
    end