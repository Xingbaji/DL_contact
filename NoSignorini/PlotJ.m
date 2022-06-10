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
axis on;  %���������Ὺ��
xlabel('u_\nu')
ylabel('j ( u_\nu )')
% set(gca,'xtick',[0,0.006,0.01]);%gca�ǵ�ǰ������ľ����xtick��ʾ��Ҫ����x��̶�Ҫ��ʾ��λ��
% set(gca,'xticklabel',{0,'u_\nu^1',' g'}); %xticklabel��ʾ���ÿ̶�����ʾ�Ķ���������Ϊϣ����ʾ��ʵ��ֵ
% set(gca,'ytick',[]);%����ط�ע��y���Ǵ����������ģ�0��������
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
    