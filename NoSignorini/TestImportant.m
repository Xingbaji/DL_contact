%��Ҫ��test���򣺻�ͼ������ÿ��С�����������е�ֵ����ɫ��ʾ����
TRI = MyMeshRectangleLeftD1(0,1,0,1,2,2);
sigma=[1,2,3,4,5,6,7,8];
big=max(sigma);
Nt=size(TRI.Elements,1);

for i=1:Nt
    coords=MyGetNodes(TRI,i);
    mesh(coords(:,1),coords(:,2),sigma(i))  
    %fill(coords(:,1),coords(:,2),[sigma(i)/big,sigma(i)/big,sigma(i)/big])  %ֵС����ɫ�ֵ�����ɫǳ
    hold on
    %polyin = polyshape(coords(:,1),coords(:,2));
    %TFin = isinterior(polyin,z)
end

data=xlsread('����һ���ѽ�����Ŀ��������.xls');
x=data(:,1);
y=data(:,2);
z=data(:,3);
 
S = 10; %�����Ĵ�С/�ߴ�
scatter(x,y,S,z,'filled') %filled��ʾ����ʵ�ĵ㣬ȱʡ��Ϊ���ĵ�
xlabel('γ�ȣ��㣩')
ylabel('���ȣ��㣩')
grid on
h = colorbar;
set(get(h,'label'),'string','����۸� (Ԫ)');%����ɫ������
 
xlim([22.4931 23.8784]) %����������̶�ȡֵ��Χ
ylim([112.6833 114.5130])



