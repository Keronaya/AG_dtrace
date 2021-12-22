function [pro0,pro1,pro2,pro3]=genrate_p(p,p3)
I=eye(p);
diag_1([1:p-1],1)=0.2;
diag_2([1:p-2],1)=0.2;
diag_3([1:p-3],1)=0.2;
diag_4([1:p-4],1)=0.2;
pro0=I;
pro1=I+diag(diag_1,1)+diag(diag_2,2)+diag(diag_1,-1)+diag(diag_2,-2);
pro2=I+diag(diag_1,1)+diag(diag_2,2)+diag(diag_3,3)+diag(diag_4,4)+diag(diag_1,-1)+diag(diag_2,-2)+diag(diag_3,-3)+diag(diag_4,-4);
pro3=eye(p3);
if rem(sqrt(p3),1)==0
    for i=1:p3
        if mod(i,sqrt(p3))~=0 &&i<p3
            pro3(i,i+1)=0.2;
        end
        if i+sqrt(p3)<=p3
            pro3(i,i+sqrt(p3))=0.2;
        end
    end
else
    'p3 require p3^(1/2) be an integer'
end
end

