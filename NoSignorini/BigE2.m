function M2=BigE2(sigma,coeE,coeK)
%给定网格划分T、按小三角形序号的分片常数sigma：Nt*4矩阵、BigE作用下的系数coeE与coeK
%希望能得到BigE2作用sigma后的分片常数M2：Nt*4矩阵

% 得到网格中自由顶点的数量、小三角形的数量，并初始化M

coe1=coeE*coeK/(1-coeK^2);
coe2=coeE/(1+coeK);

M2=zeros(4,1); %要得到的\BigE（\varepsion(u)）在网格的每个小三角形是分片常数的

M2(1)=coe1*( sigma(1)+sigma(4) ) +coe2*sigma(1);
M2(2)=coe2*sigma(2);
M2(3)=coe2*sigma(3);
M2(4)=coe1*( sigma(1)+sigma(4) ) +coe2*sigma(4);


