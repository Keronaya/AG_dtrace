

clear, clc
A = rand(3e3);
B = rand(3e3);
num = 3;
T = zeros(8,num);
for i = 1:num
    % Test1 乘法
    clc
    disp('^-------')
    disp(['第' sprintf('%4i',i) ' 轮乘法测试中...'])
    tic, X1 = A*B; T(1,i) = toc;
    % Test2 稀疏矩阵
    clc
    disp('-^------')
    disp(['第' sprintf('%4i',i) ' 轮稀疏矩阵测试中...'])
    tic, X2 = sparse(A); T(2,i) = toc;
    % Test3 逆矩阵
    clc
    disp('--^-----')
    disp(['第' sprintf('%4i',i) ' 轮逆矩阵测试中...'])
    tic, X3 = inv(A); T(3,i) = toc;
    % Test4 快速傅里叶
    clc
    disp('---^----')
    disp(['第' sprintf('%4i',i) ' 轮快速傅里叶测试中...'])
    tic, X4 = fft(A); T(4,i) = toc;
    % Test5 LU分解
    clc
    disp('----^---')
    disp(['第' sprintf('%4i',i) ' 轮LU分解测试中...'])
    tic, [L5,U5,P5] = lu(A); T(5,i) = toc;
    % Test6 QR分解
    clc
    disp('-----^--')
    disp(['第' sprintf('%4i',i) ' 轮QR分解测试中...'])
    tic, X6 = qr(A); T(6,i) = toc;
    % Test7 奇异值分解
    clc
    disp('------^-')
    disp(['第' sprintf('%4i',i) ' 轮奇异值分解测试中...'])
    tic, [U7,S7,V7] = svd(A); T(7,i) = toc;
    % Test8 特征值与特征向量
    clc
    disp('-------^')
    disp(['第' sprintf('%4i',i) ' 轮特征值与特征向量测试中...'])
    tic, [V8,D8] = eig(A); T(8,i) = toc;
end
clc
% 各项测试平均时间
t = sum(T,2)./num;
disp(['Multiplication :   ' sprintf('%6f',t(1))])
disp(['Sparse         :   ' sprintf('%6f',t(2))])
disp(['Inverse        :   ' sprintf('%6f',t(3))])
disp(['FFT            :   ' sprintf('%6f',t(4))])
disp(['LU             :   ' sprintf('%6f',t(5))])
disp(['QR             :   ' sprintf('%6f',t(6))])
disp(['SVD            :   ' sprintf('%6f',t(7))])
disp(['Eigen          :   ' sprintf('%6f',t(8))])
% Total
total = sum(t);
disp('----------------------------')
disp(['Total          :   ' sprintf('%6f',total)])