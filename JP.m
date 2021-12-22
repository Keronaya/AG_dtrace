function [U T]=JP(Sigma,Theta,p)
I=eye(p);
T=zeros(p);
Gamma=(1/2)*(Sigma*Theta+Theta*Sigma);
Gamma_1=Gamma-I;
A=Gamma_1;
T=Theta-A;                                    
U=zeros(p);
for i=1:p
    for j=1:p
        U(i,j)=(1/4)*Sigma(i,i)*Theta(j,j)+(1/4)*Sigma(j,j)*Theta(i,i)+(1/4)*Sigma(i,j)*Theta(i,j)+(1/2)*I(i,j)+1/2; %%% \sigma
    end
end
end
 %%%%% 此函数为debiased 估计量