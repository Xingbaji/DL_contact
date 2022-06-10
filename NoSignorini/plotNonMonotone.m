x=-0.005 : 0.001 : 0.01;
N=length(x);
y=zeros(N,1);

for i=1:N
    if x(i)<=0
        y(i)=0;
    elseif x(i)<=0.006
        y(i)=300*x(i);
    else
        y(i)=300*0.006-100*(x(i)-0.006);
    end
end

figure1 = figure;
%annotation1 = annotation(figure1,'arrow',[0.131 0.131],[0.70 0.85]);
annotation1 = annotation(figure1,'arrow',[0.325 0.325],[0.11 0.85]);%y轴箭头
annotation2 = annotation(figure1,'arrow',[0.131 0.85],[0.11 0.11]); %x轴箭头


plot(x,y,'r','linewidth',1.5)
axis on;  %设置坐标轴开启
set(gca,'xtick',[0,0.006,0.01]);%gca是当前坐标轴的句柄，xtick表示我要设置x轴刻度要显示的位置
set(gca,'xticklabel',{0,'u_\nu^1',' g'}); %xticklabel表示设置刻度上显示的东西，后面为希望显示的实际值
set(gca,'ytick',[]);%这个地方注意y轴是从上往下数的，0在最上面
set(gca,'yticklabel',{});

text(0.013,-0.15,'u_\nu')
text(-0.001 , 2.2 , '-\sigma_\nu')

axis([-0.005 0.015 0 2.7])
hold on
xx=[0.01,0.01];
yy=[y(N),2.7];
plot(xx,yy,'r','linewidth',1.5)
    