function [xs] = algo1(p,lambda,T,U )%uΪЭ���� %t��ƫ����
xs=zeros(p);
for i=1:p
    for j=1:p
        if (abs(T(i,j))>lambda*sqrt(U(i,j)))
            xs(i,j)=T(i,j);
        else
            xs(i,j)=0;
        end
    end
end
end

