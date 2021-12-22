function [pre]=Examples(p)% Sigman 样本协方差，Sigma0 样本方差
pre=zeros(p);
for i=1:p
    pre(i,i)=1;
end
for i=1:p-1
    pre(i,i+1)=0.4;
end
for i=2:p
    pre(i,i-1)=0.4;
end
end
