

clear, clc
A = rand(3e3);
B = rand(3e3);
num = 3;
T = zeros(8,num);
for i = 1:num
    % Test1 �˷�
    clc
    disp('^-------')
    disp(['��' sprintf('%4i',i) ' �ֳ˷�������...'])
    tic, X1 = A*B; T(1,i) = toc;
    % Test2 ϡ�����
    clc
    disp('-^------')
    disp(['��' sprintf('%4i',i) ' ��ϡ����������...'])
    tic, X2 = sparse(A); T(2,i) = toc;
    % Test3 �����
    clc
    disp('--^-----')
    disp(['��' sprintf('%4i',i) ' ������������...'])
    tic, X3 = inv(A); T(3,i) = toc;
    % Test4 ���ٸ���Ҷ
    clc
    disp('---^----')
    disp(['��' sprintf('%4i',i) ' �ֿ��ٸ���Ҷ������...'])
    tic, X4 = fft(A); T(4,i) = toc;
    % Test5 LU�ֽ�
    clc
    disp('----^---')
    disp(['��' sprintf('%4i',i) ' ��LU�ֽ������...'])
    tic, [L5,U5,P5] = lu(A); T(5,i) = toc;
    % Test6 QR�ֽ�
    clc
    disp('-----^--')
    disp(['��' sprintf('%4i',i) ' ��QR�ֽ������...'])
    tic, X6 = qr(A); T(6,i) = toc;
    % Test7 ����ֵ�ֽ�
    clc
    disp('------^-')
    disp(['��' sprintf('%4i',i) ' ������ֵ�ֽ������...'])
    tic, [U7,S7,V7] = svd(A); T(7,i) = toc;
    % Test8 ����ֵ����������
    clc
    disp('-------^')
    disp(['��' sprintf('%4i',i) ' ������ֵ����������������...'])
    tic, [V8,D8] = eig(A); T(8,i) = toc;
end
clc
% �������ƽ��ʱ��
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