%此程序是为了得到 \|u_n\|_V 和\|\sigma\|_Q的值

X = [0,1];                         
Y = [0,1]; 
Time=1; 

intervals2=64;total_time_step2=64;  
T = MyMeshRectangleLeftD1(X(1),X(2),Y(1),Y(2),intervals2,intervals2);
FNd2=intervals2*(intervals2+1);

 %密网格下作为标准解
load('u64_64.mat')     
load('sigma64_64.mat')
Uref=u64_64;
Sigmaref=sigma64_64;


NormSigma=0.5*(1/intervals2)^2 * sum( sqrt(Sigmaref(:,1,total_time_step2+1).^2 +...
    Sigmaref(:,2,total_time_step2+1).^2 +Sigmaref(:,3,total_time_step2+1).^2 +Sigmaref(:,4,total_time_step2+1).^2) );
NormU=MyErrorStiffness2(T,Uref(1:2*FNd2,total_time_step2+1));


MaxSigma=zeros(total_time_step2+1,1);
MaxU=zeros(total_time_step2+1,1);
MaxPlus=zeros(total_time_step2+1,1);
for stepp=1:total_time_step2+1
    SigmaCha=Sigmaref(:,:,stepp); 
    MaxSigma(stepp)=0.5*(1/intervals2)^2 * sum( sqrt(SigmaCha(:,1).^2 +SigmaCha(:,2).^2 +SigmaCha(:,3).^2 +SigmaCha(:,4).^2) );   %一般的sigma范数
    %U的误差
    MaxU(stepp)=MyErrorStiffness2(T,Uref(1:2*FNd2,total_time_step2));
    
    MaxPlus(stepp)=MaxU(stepp)+MaxSigma(stepp);
end
    