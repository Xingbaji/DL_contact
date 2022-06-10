% �����棬��ȷ��

X = [0,1];                         
Y = [0,1]; 
Time=1; 

intervals1=4;  %������1
intervals2=64; %������2
Rate=round(intervals2/intervals1);
iteration=log2(Rate); %ÿ�ν�������2�����ܣ���Ҫ����iteration��

%�������µ���ֵ��
[U,~]=MyContactVStatic(intervals1);   
% load('u32.mat')     
% U=u32;

%����������Ϊ��׼��
%[Uref,~]=MyContactVStatic(intervals2);   
load('u64.mat')     
Uref=u64;



%for stepp=1:total_time_step1+1

%���������¶�U���ֽ��в�ֵ���õ�UInterp
FNd1=intervals1*(intervals1+1);
FNd2=intervals2*(intervals2+1);
UX=U(1:FNd1);UY=U(FNd1+1:2*FNd1);
UInterpX=zeros(FNd2,1);
UInterpY=zeros(FNd2,1);


for kx=1:(intervals1+1)
    for jx=1:intervals1
        if jx==1
            for ix=1:Rate
                UInterpX( (kx-1)*intervals2*Rate+ix )=ix/Rate*UX( (kx-1)*intervals1+jx);
                UInterpY( (kx-1)*intervals2*Rate+ix )=ix/Rate*UY( (kx-1)*intervals1+jx);
            end
        else
            for ix=1:Rate
                UInterpX( (kx-1)*intervals2*Rate+Rate*(jx-1)+ix ) =(1-ix/Rate) *UX( (kx-1)*intervals1+jx-1 ) + ix/Rate *UX( (kx-1)*intervals1+jx );
                UInterpY( (kx-1)*intervals2*Rate+Rate*(jx-1)+ix ) =(1-ix/Rate) *UY( (kx-1)*intervals1+jx-1 ) + ix/Rate *UY( (kx-1)*intervals1+jx );
            end
        end
    end 
end

for kx=1:intervals1
    for jx=1:(Rate-1)
        for ix=1:intervals2
            UInterpX( (kx-1)*intervals2*Rate+ jx*intervals2 +ix  )=(1-jx/Rate) *UInterpX( (kx-1)*intervals2*Rate+ix ) +...
                jx/Rate *UInterpX( kx*intervals2*Rate+ix );
            UInterpY( (kx-1)*intervals2*Rate+ jx*intervals2 +ix  )=(1-jx/Rate) *UInterpY( (kx-1)*intervals2*Rate+ix ) +...
                jx/Rate *UInterpY( kx*intervals2*Rate+ix );
        end
    end
end
UInterp=[UInterpX;UInterpY];

T = MyMeshRectangleLeftD1(X(1),X(2),Y(1),Y(2),intervals2,intervals2);


%U�����
%ErrorU=MyErrorStiffness2(T,UInterp-Uref(1:2*FNd2));
ErrorU=sqrt( (UInterp-Uref(1:2*FNd2))' * MyOriginalStiff2(T) ...   %������������ע������ġ���������������������
    * (UInterp-Uref(1:2*FNd2)) / (2*intervals2^2));


%BiaoZhunU=MyErrorStiffness2(T,Uref(1:2*FNd2));
BiaoZhunU=sqrt( Uref(1:2*FNd2)' * MyOriginalStiff2(T) ...   %������������ע������ġ���������������������
    * Uref(1:2*FNd2) / (2*intervals2^2));
