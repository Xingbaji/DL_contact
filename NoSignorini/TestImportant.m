%重要的test程序：画图——将每个小三角形网格中的值用颜色表示出来
TRI = MyMeshRectangleLeftD1(0,1,0,1,2,2);
sigma=[1,2,3,4,5,6,7,8];
big=max(sigma);
Nt=size(TRI.Elements,1);

for i=1:Nt
    coords=MyGetNodes(TRI,i);
    mesh(coords(:,1),coords(:,2),sigma(i))  
    %fill(coords(:,1),coords(:,2),[sigma(i)/big,sigma(i)/big,sigma(i)/big])  %值小的颜色深，值大的颜色浅
    hold on
    %polyin = polyshape(coords(:,1),coords(:,2));
    %TFin = isinterior(polyin,z)
end

data=xlsread('附件一：已结束项目任务数据.xls');
x=data(:,1);
y=data(:,2);
z=data(:,3);
 
S = 10; %坐标点的大小/尺寸
scatter(x,y,S,z,'filled') %filled表示点是实心点，缺省则为空心点
xlabel('纬度（°）')
ylabel('经度（°）')
grid on
h = colorbar;
set(get(h,'label'),'string','任务价格 (元)');%给颜色栏命名
 
xlim([22.4931 23.8784]) %设置坐标轴刻度取值范围
ylim([112.6833 114.5130])



