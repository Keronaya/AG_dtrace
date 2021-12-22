function [s0,it,times]=servers_do_SP(s0,rho,l_n,MAXIT,TOL1,eps,p,N,times)
%sigam_hat   1
%theta       2
%theta0      3
%theta1      4
%´ólambda0   5
%´ólambda1   6
%theta 6
EPS=eps*eye(p);
it=0;

MAXIT_in=1000;
TOL = 10^(-7);
for iter=1:MAXIT
    if iter==1
        for i=1:N
            [s0{i,2},times(iter,i)]=argmins(s0{i,1},rho,l_n,MAXIT_in,TOL,EPS,p);
            [U,T]=JP(s0{i,1},s0{i,2},p);
            s0{i,4}= algo1(p,l_n,T,U );
            %s0{i,2}=s0{i,4};
        end
        iter=iter+1;
    end
    for i=1:2:N
        if i==1
            %s0{i,2}=G_1(s0{i,1}+2*rho*I,I+rho*s0{i+1,2}-s0{i,3});%theta
            s0{i,4}=(s0{i,4}+s0{i+1,4})/2;Thetaj=(s0{i,2}+s0{i+1,2})/2;Sigmaj=(s0{i,1}+s0{i+1,1})/2;
            s0{i,2}=algo2(p,l_n/2,s0{i,4},Thetaj,Sigmaj); 
            s0{i,4}=s0{i,2};
            %s0{i,3}=s0{i,3}+rho*(s0{i,2}-s0{i+1,2});%lambda1
        else
            
            s0{i,4}=(s0{i,4}+s0{i-1,4}+s0{i+1,4})/3;Thetaj=(s0{i,2}+s0{i-1,2}+s0{i+1,2})/3;Sigmaj=(s0{i,1}+s0{i-1,1}+s0{i+1,1})/3;
            s0{i,2}=algo2(p,l_n/3,s0{i,4},Thetaj,Sigmaj);
            s0{i,4}=s0{i,2};
            %s0{i,3}=s0{i,3}+rho*(s0{i,2}-s0{i+1,2});%lambda1
        end
    end
    
    
    %tail update
    for i=2:2:N
        if i==N
            s0{i,4}=(s0{i,4}+s0{i-1,4})/2;Thetaj=(s0{i,2}+s0{i-1,2})/2;Sigmaj=(s0{i,1}+s0{i-1,1})/2;
            s0{i,2}=algo2(p,l_n/2,s0{i,4},Thetaj,Sigmaj);
            s0{i,4}=s0{i,2};
            %s0{i,3}=s0{i,3}+rho*(s0{i,2}-s0{i+1,2});%lambda1
        else
            s0{i,4}=(s0{i,4}+s0{i-1,4}+s0{i+1,4})/3;Thetaj=(s0{i,2}+s0{i-1,2}+s0{i+1,2})/3;Sigmaj=(s0{i,1}+s0{i-1,1}+s0{i+1,1})/3;
            s0{i,2}=algo2(p,l_n/3,s0{i,4},Thetaj,Sigmaj);
            s0{i,4}=s0{i,2};
            %s0{i,3}=s0{i,3}+rho*(s0{i,2}-s0{i+1,2});%lambda1
        end %(s_hat,theta0,theta1,lambda0,lambda1,rho,l_n,MAXIT_in,TOL,EPS,p)
    end
    
    %update head
    for i=1:2:N
        if i==1
            %s0{i,2}=G_1(s0{i,1}+2*rho*I,I+rho*s0{i+1,2}-s0{i,3});%theta
            [s0{i,2},s0{i,3},times(iter,i)]=argmins1_np(s0{i,1},s0{i,2},s0{i+1,2},s0{i,3},rho,l_n,MAXIT_in,TOL,EPS,p);
            s0{i,4}=(s0{i,4}+s0{i+1,4})/2;Thetaj=(s0{i,2}+s0{i+1,2})/2;Sigmaj=(s0{i,1}+s0{i+1,1})/2;
            s0{i,2}=algo2(p,l_n/2,s0{i,4},Thetaj,Sigmaj); 
            s0{i,4}=s0{i,2};
            %s0{i,3}=s0{i,3}+rho*(s0{i,2}-s0{i+1,2});%lambda1
        else
            [s0{i,2},s0{i-1,3},s0{i,3},times(iter,i)]=argmins_np(s0{i,1},s0{i,2},s0{i-1,2},s0{i+1,2},s0{i-1,3},s0{i,3},rho,l_n,MAXIT_in,TOL,EPS,p);
            s0{i,4}=(s0{i,4}+s0{i-1,4}+s0{i+1,4})/3;Thetaj=(s0{i,2}+s0{i-1,2}+s0{i+1,2})/3;Sigmaj=(s0{i,1}+s0{i-1,1}+s0{i+1,1})/3;
            s0{i,2}=algo2(p,l_n/3,s0{i,4},Thetaj,Sigmaj);
            s0{i,4}=s0{i,2};
            %s0{i,3}=s0{i,3}+rho*(s0{i,2}-s0{i+1,2});%lambda1
        end
    end
    juedgeif=[norm(s0{N-1,2}-s0{N,2},'fro')/max([1,norm(s0{N,2},'fro'),norm(s0{N-1,2},'fro')]),
        norm(s0{N-1,2}-s0{N-2,2},'fro')/max([1,norm(s0{N-1,2},'fro'),norm(s0{N-2,2},'fro')]),
        norm(s0{N-2,2}-s0{N-3,2},'fro')/max([1,norm(s0{N-2,2},'fro'),norm(s0{N-3,2},'fro')])]
    it=iter
    if juedgeif<TOL1 
        'success'
        return
    end
    %'head update'
    %'µü´ú'
    
    
    %tail update
    for i=2:2:N
        if i==N
            [s0{i,2},s0{i-1,3},times(iter,i)]=argminsn_np(s0{i,1},s0{i,2},s0{i-1,2},s0{i-1,3},rho,l_n,MAXIT_in,TOL,EPS,p);
            s0{i,4}=(s0{i,4}+s0{i-1,4})/2;Thetaj=(s0{i,2}+s0{i-1,2})/2;Sigmaj=(s0{i,1}+s0{i-1,1})/2;
            s0{i,2}=algo2(p,l_n/2,s0{i,4},Thetaj,Sigmaj);
            s0{i,4}=s0{i,2};
            %s0{i,3}=s0{i,3}+rho*(s0{i,2}-s0{i+1,2});%lambda1
        else
            [s0{i,2},s0{i-1,3},s0{i,3},times(iter,i)]=argmins_np(s0{i,1},s0{i,2},s0{i-1,2},s0{i+1,2},s0{i-1,3},s0{i,3},rho,l_n,MAXIT_in,TOL,EPS,p);
            s0{i,4}=(s0{i,4}+s0{i-1,4}+s0{i+1,4})/3;Thetaj=(s0{i,2}+s0{i-1,2}+s0{i+1,2})/3;Sigmaj=(s0{i,1}+s0{i-1,1}+s0{i+1,1})/3;
            s0{i,2}=algo2(p,l_n/3,s0{i,4},Thetaj,Sigmaj);
            s0{i,4}=s0{i,2};
            %s0{i,3}=s0{i,3}+rho*(s0{i,2}-s0{i+1,2});%lambda1
        end %(s_hat,theta0,theta1,lambda0,lambda1,rho,l_n,MAXIT_in,TOL,EPS,p)
    end
    %'tail update'   
end
end






