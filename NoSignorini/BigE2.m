function M2=BigE2(sigma,coeE,coeK)
%�������񻮷�T����С��������ŵķ�Ƭ����sigma��Nt*4����BigE�����µ�ϵ��coeE��coeK
%ϣ���ܵõ�BigE2����sigma��ķ�Ƭ����M2��Nt*4����

% �õ����������ɶ����������С�����ε�����������ʼ��M

coe1=coeE*coeK/(1-coeK^2);
coe2=coeE/(1+coeK);

M2=zeros(4,1); %Ҫ�õ���\BigE��\varepsion(u)���������ÿ��С�������Ƿ�Ƭ������

M2(1)=coe1*( sigma(1)+sigma(4) ) +coe2*sigma(1);
M2(2)=coe2*sigma(2);
M2(3)=coe2*sigma(3);
M2(4)=coe1*( sigma(1)+sigma(4) ) +coe2*sigma(4);


