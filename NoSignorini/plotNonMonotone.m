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
annotation1 = annotation(figure1,'arrow',[0.325 0.325],[0.11 0.85]);%y���ͷ
annotation2 = annotation(figure1,'arrow',[0.131 0.85],[0.11 0.11]); %x���ͷ


plot(x,y,'r','linewidth',1.5)
axis on;  %���������Ὺ��
set(gca,'xtick',[0,0.006,0.01]);%gca�ǵ�ǰ������ľ����xtick��ʾ��Ҫ����x��̶�Ҫ��ʾ��λ��
set(gca,'xticklabel',{0,'u_\nu^1',' g'}); %xticklabel��ʾ���ÿ̶�����ʾ�Ķ���������Ϊϣ����ʾ��ʵ��ֵ
set(gca,'ytick',[]);%����ط�ע��y���Ǵ����������ģ�0��������
set(gca,'yticklabel',{});

text(0.013,-0.15,'u_\nu')
text(-0.001 , 2.2 , '-\sigma_\nu')

axis([-0.005 0.015 0 2.7])
hold on
xx=[0.01,0.01];
yy=[y(N),2.7];
plot(xx,yy,'r','linewidth',1.5)
    