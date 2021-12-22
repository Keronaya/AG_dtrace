function [s0,l_n,itout,itin,times_in]=lla(s0,rho,l_n,MAXIT,TOL,EPS,p12,N,times)
a=3;
TOLsp=10^(-6);
for i=1:10000
    th_pre=s0{1,2};
    [s0,itin,times_in]=servers_do_SP(s0,rho,l_n,MAXIT,TOLsp,EPS,p12,N,times);
    l_n=update_ln(l_n,abs(s0{1,2}),a);
    if norm(s0{1,2}-th_pre,Inf)<=TOL
        itout=i;
        return
    end
end
    
end