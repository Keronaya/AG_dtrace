function [th_n_kp1]=G_1(A,B)
D=size(A,1);
[U_A,S_A]=eig(A);
U_A=real(U_A);
S_A=real(S_A);

%S_A = fliplr(S_A);
%V_A=real(V_A);
temp1=repmat(diag(S_A)./2,1,D);
temp1=temp1+temp1';
th_n_kp1=U_A*((U_A'*B*U_A)./temp1)*U_A';


