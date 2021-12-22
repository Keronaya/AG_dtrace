function [th1_n_kp1]=P_1(X,eps)
%D=size(X,2);
[U_X,S_X]=eig(X);
U_X=real(U_X);
S_X=real(S_X);
th1_n_kp1=U_X*max(S_X,eps)*U_X';