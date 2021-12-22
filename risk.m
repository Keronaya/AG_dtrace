function [fro,op,mat,tp,tn,FPR,FNR]=risk(theta,pro)
fro=norm(theta-pro,'fro');
op=norm(theta-pro,2);
mat=norm(theta-pro,inf);
theta=round(theta,6);
p=size(pro,1);
tp=sum(sum(logical(theta)&logical(pro)))/max(sum(sum(logical(pro)==1)),sum(sum(logical(theta)==1)));
tn=sum(sum(logical(theta==0)&logical(pro==0)))/max(sum(sum(logical(pro)==0)),sum(sum(logical(theta)==0)));
FPR=sum(sum(logical(theta~=0)&logical(pro==0)))/sum(sum(logical(pro==0)));
FNR=sum(sum(logical(theta==0)&logical(pro~=0)))/sum(sum(logical(pro~=0)));
end