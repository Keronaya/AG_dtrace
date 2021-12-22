function [theta,itin]=argmins(s_hat,rho,l_n,MAXIT_in,TOL,EPS,p)
I=eye(p);
theta=inv(diag(diag(s_hat)));
theta0=theta;
theta1=theta0;
lambda0=zeros(p);
lambda1=zeros(p);
itin=0;
for iter=1:MAXIT_in
    theta_o=theta;
    
    theta0_o=theta0;
    theta1_o=theta1;
    %lambda0_o=lambda0;
    %lambda1_o=lambda1;
    %theta
    theta=G_1(s_hat+2*rho*I,I+rho*(theta0+theta1)-lambda0-lambda1);
    theta0=S_1(theta+lambda0/rho,l_n/rho);%theta0
    theta1=P_1(theta+lambda1/rho,EPS);%theta1
    lambda0=lambda0+rho*(theta-theta0); %lambda0    may be it can put on the other side
    lambda1=lambda1+rho*(theta-theta1); %lambda1
    
    if norm(theta-theta_o,'fro')/max([1,norm(theta_o,'fro'),norm(theta,'fro')])<TOL && norm(theta0-theta0_o,'fro')/max([1,norm(theta0_o,'fro'),norm(theta0,'fro')])<TOL && norm(theta1-theta1_o,'fro')/max([1,norm(theta1_o,'fro'),norm(theta1,'fro')])<TOL
        itin=iter;
        return
    end
   % 'argmins'
    %iter
end


end
