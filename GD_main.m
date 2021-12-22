clc;
% clear all;
N=4;%servers numbers
sample=400;
All_server_sample=N*sample;
%l_n=0.1; %Сlambda
MAXIT = 1000;
TOL = 10^(-6);
EPS = 0;
p12 = 500 ;%dimension
p3 = 484;
rho=1;
l_n=sqrt(log(p12)/sample);
l_n_3=sqrt(log(p3)/sample);
time={};
risks0=zeros(100,5);
risks1=zeros(100,5);
risks2=zeros(100,5);
risks3=zeros(100,5);
itall=zeros(100,4);
% X0= csvread('C:\Users\16117\Desktop\codeEx\Python\GADMM\datax0.csv');
%X1= csvread('C:\Users\16117\Desktop\codeEx\Python\GADMM\datax1.csv');
% X2= csvread('C:\Users\16117\Desktop\codeEx\Python\GADMM\datax2.csv');
%X3_b= csvread('C:\Users\16117\Desktop\codeEx\Python\GADMM\datax3.csv');
% pr0= csvread('C:\Users\16117\Desktop\codeEx\Python\GADMM\pr0.csv');
%pr1= csvread('C:\Users\16117\Desktop\codeEx\Python\GADMM\pr1.csv');
%pr2= csvread('C:\Users\16117\Desktop\codeEx\Python\GADMM\pr2.csv');
%pr3_b= csvread('C:\Users\16117\Desktop\codeEx\Python\GADMM\pr3.csv');
%[pro0,pro1,pro2,pro3]=genrate_p(p12,p3);
u=zeros(p12,1)';
p_w=Examples(p12);
R=inv(p_w);
%x3=multivrandn(u,pro3,sample);
%pro3=(pro3+pro3')/2;
%server=servers_data(X0,sample,N);
%[x0,x1,x2,x3]=genrate_x(pro0,pro1,pro2,pro3,All_server_sample);
% server1=servers_data(x1,sample,N);
% server2=servers_data(X2,sample,N);
for i=1:5
    %[x0,x1,x2,x3]=genrate_x(pro0,pro1,pro2,pro3,All_server_sample);
    %server0=servers_data(x0,sample,N);
     x_w=multivrandn(u,R,All_server_sample);
     server1=servers_data(x_w,sample,N);
%     server2=servers_data(x2,sample,N);
%     server3=servers_data(x3,sample,N);
    %times0=zeros(MAXIT,N);
     times1=zeros(MAXIT,N);
%     times2=zeros(MAXIT,N);
%     times3=zeros(MAXIT,N);
    %[aim0,it0,times0]=servers_do_SP(server0,rho,l_n+0.01*i,MAXIT,TOL,EPS,p12,N,times0);
    [aim1,it1,times1]=servers_do_SP(server1,rho,l_n,MAXIT,TOL,EPS,p12,N,times1);
    %[aim2,it2,times2]=servers_do_SP(server2,rho,l_n,MAXIT,TOL,EPS,p12,N,times2);
    %[aim3,it3,times3]=servers_do_SP(server3,rho,l_n_3,MAXIT,TOL,EPS,p3,N,times3);
    %time{i,1}=times0;
     time{i,2}=times1;
%     time{i,3}=times2;
%     time{i,4}=times3;
    %itall(i,1)=it0;
     itall(i,2)=it1;
%     itall(i,3)=it2;
%     itall(i,4)=it3;
    %[risks0(i,1),risks0(i,2),risks0(i,3),risks0(i,4),risks0(i,5)]=risk(aim0{1,2},pro0);
    [risks1(i,1),risks1(i,2),risks1(i,3),risks1(i,4),risks1(i,5)]=risk(aim1{3,2},p_w);
    %[risks2(i,1),risks2(i,2),risks2(i,3),risks2(i,4),risks2(i,5)]=risk(aim2{1,2},pro2);
    %[risks3(i,1),risks3(i,2),risks3(i,3),risks3(i,4),risks3(i,5)]=risk(aim3{1,2},pro3);
end

% tic
% [aim1,it1,times]=servers_do_SP(server1,rho,0.1269,MAXIT,TOL,EPS,p12,N,times);
% toc
% tic
% [target1,l_n1,itout1,itin1,times_in1]=lla(server1,rho,0.078,MAXIT,TOL,EPS,p12,N,times);
% toc


% [thetaa,theta0a,ita]=alltest(X1,rho,0.13,MAXIT,TOL,EPS,p12);