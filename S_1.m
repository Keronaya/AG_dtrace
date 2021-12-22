function [th0_n_kp1]=S_1(A,l_n_r)
%I=eye(size(A,1));
Aii=diag(diag(A)); %主对角备用
Aij=(A-l_n_r).*logical(A>l_n_r)+(A+l_n_r).*logical(A<-l_n_r);
Aij=Aij-diag(diag(Aij));
th0_n_kp1=Aii+Aij;


