function [x0,x1,x2,x3]=genrate_x(pro0,pro1,pro2,pro3,allsample)
p012=size(pro0,1);
p3=size(pro3,1);
sigma0=inv(pro0);
sigma1=inv(pro1);
sigma2=inv(pro2);
sigma3=inv(pro3);
x0=mvnrnd(zeros(1,p012),sigma0,allsample);
x1=mvnrnd(zeros(1,p012),sigma1,allsample);
x2=mvnrnd(zeros(1,p012),sigma2,allsample);
x3=mvnrnd(zeros(1,p3),sigma3,allsample);
end