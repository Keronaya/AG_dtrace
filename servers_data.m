function [server]=servers_data(X,sample,N)
p=size(X,2);
%server={zeros(p),zeros(p),zeros(p),zeros(p);zeros(p),zeros(p),zeros(p);zeros(p),zeros(p),zeros(p);zeros(p),zeros(p),zeros(p);};
for i=1:N      %服务器参数初始化
    server{i,1}=cov(X((i-1)*sample+1:i*sample,:));%X数据分布在各个服务器
    server{i,2}=zeros(p);%theta
    server{i,3}=zeros(p);%lambda
    server{i,4}=zeros(p);
end
end
