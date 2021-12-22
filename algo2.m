function [xsf] =algo2(p,tau,xsj,Thetaj,Sigmaj)
xsf=zeros(p);U2=zeros(p);I=eye(p);
for i=1:p
    for j=1:p
     U2(i,j)=(1/4)*Sigmaj(i,i)*Thetaj(j,j)+(1/4)*Sigmaj(j,j)*Thetaj(i,i)+(1/4)*Sigmaj(i,j)*Thetaj(i,j)+(1/2)*I(i,j)+1/2;
    end
end 
for i=1:p
    for j=1:p
        if (abs(xsj(i,j))>tau*sqrt(U2(i,j)))
            xsf(i,j)=xsj(i,j);
        else
            xsf(i,j)=0;
        end
    end
end
end

