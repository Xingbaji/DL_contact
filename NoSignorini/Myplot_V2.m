% 用h=1/64,k(tau)=1/128作为标准参考解；
%对比h=1/2, 1/4, 1/8, 1/16, 1/32
b=[8.9609, 5.4782, 3.0784, 1.5060, 0.5657 ];
a=[1,1/2,1/4,1/8,1/16];
plot(a, b,'-*')
%title('Error ~ (h+\tau)')
xlabel('h+\tau')
ylabel('Error')
