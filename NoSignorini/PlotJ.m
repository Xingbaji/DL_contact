x=-0.002 : 0.0005 : 0.012;
N=length(x);
y=zeros(N,1);

for i=1:N
    if x(i)<=0
        y(i)=0;
    elseif x(i)<=0.006
        y(i)=150*x(i)^2;
    elseif x(i) <=0.01
        y(i)=-50*x(i)^2+2.4*x(i)-0.0072;
    else
        y(i)=150*x(i)^2 - 1.6*x(i)+0.0128;
    end
end


plot(x,y,'r','linewidth',1.5)
axis on;  %设置坐标轴开启
xlabel('u_\nu')
ylabel('j ( u_\nu )')
% set(gca,'xtick',[0,0.006,0.01]);%gca是当前坐标轴的句柄，xtick表示我要设置x轴刻度要显示的位置
% set(gca,'xticklabel',{0,'u_\nu^1',' g'}); %xticklabel表示设置刻度上显示的东西，后面为希望显示的实际值
% set(gca,'ytick',[]);%这个地方注意y轴是从上往下数的，0在最上面
% set(gca,'yticklabel',{});
% 
% text(0.013,-0.15,'u_\nu')
% text(-0.001 , 2.2 , '-\sigma_\nu')
% 
% axis([-0.005 0.015 0 2.7])
% hold on
% xx=[0.01,0.01];
% yy=[y(N),2.7];
% plot(xx,yy,'r','linewidth',1.5)
    