% 静态弹性接触问题
%无Sigonorini条件


function [solution_u,iteration,state]=MyContactVStatic(intervals)

%弹性方程数据
X = [0,1];
Y = [0,1];
coeE=1000;
coeK=0.4;
FExternal=-8;
%FExternal=-16; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%coe1=coeE*coeK/(1-coeK^2);
%coe2=coeE/(1+coeK);

%网格数据    
nx=intervals;
ny=intervals;
spaceh = 1/intervals;
TRI = MyMeshRectangleLeftD1(X(1),X(2),Y(1),Y(2),nx,ny);  %用到子函数《MyMeshRectangleLeftD1》
FNd=length(TRI.FNodePtrs);
Nt=size(TRI.Elements,1);

%边界Gamma3上法向次微分的数据
rrr1=0.006;rrr2=0.01; 
sharp1=300;sharp2=-100; sharp3=300;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
level1=sharp1*rrr1;
level2=sharp1*rrr1+sharp2*(rrr2-rrr1);

%求解线性方程组的格式
% The global  coefficient matrix is given by
%           K   K1   K2      K   K1
% K_whole=  T1  T2   T3    =   N          
%           
%  F_whole= F
%           FN
K=MyStiffness2(TRI,coeE,coeK);      %用到子函数《MyStiffness2》

%分配空间
solution_u=zeros(2*FNd+nx,1);
state= 1 *  ones( intervals,1 );


%边界Gamma2上受力对Load的影响
F1=zeros(2*FNd,1);
F1(2*FNd-intervals+1:end)=FExternal*spaceh;
F1(2*FNd)=F1(2*FNd)/2;

%for time_step = 2 : total_time_step+1
    complete=0;
    iteration=0;
    
    K1=sparse(zeros(2*FNd,nx));
    for i_x=1:(nx-1)
        K1( FNd+i_x, i_x :i_x+1 )=[-0.5;-0.5]*spaceh;
    end
    K1(FNd+nx,nx)=-0.5*spaceh; 
    
    %active set method loop     %%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%    %%%%%%%%%%%%%%%%%%%%
    while(complete==0 && iteration<50)   
        N=sparse(zeros(nx,2*FNd+nx));
        FN=zeros(nx,1);     
        for i_x=1:intervals
            local_1=i_x+FNd;
            local_2=i_x+2*FNd;
            if state(i_x)==1
                if i_x>1
                    N(i_x, [i_x-1+FNd, i_x+FNd, local_2] )=[0.5*sharp3,0.5*sharp3,1];                    
                else
                    N(i_x,[i_x+FNd,local_2])=[0.5*sharp3,1];
                end
                FN(i_x)=level2-sharp3*rrr2;
            elseif state(i_x)==2 
                if i_x>1
                    N(i_x, [i_x-1+FNd, i_x+FNd, local_2] )=[0.5*sharp2,0.5*sharp2,1];                     
                else
                    N(i_x,[i_x+FNd,local_2])=[0.5*sharp2,1];
                end
                FN(i_x)=level1-sharp2*rrr1;
            elseif state(i_x)==3
                if i_x>1
                    N(i_x,[i_x-1+FNd,i_x+FNd,local_2])=[0.5*sharp1,0.5*sharp1,1];                    
                else
                    N(i_x,[i_x+FNd,local_2])=[0.5*sharp1,1];
                end
                FN(i_x)=0;
            elseif state(i_x)==4 
                N(i_x,local_2)=1;
                FN(i_x)=0;
            end
        end
        
        F2 = F1; %《《《《《《《《《《《《《《《《《《《《《《《《《《《《
        K_whole=[K,K1;N];
        F_whole=[F2;FN];
        solution_u(:)=K_whole\F_whole;
        change=0;
        
        % 4 cases in total
        for i_x=1:intervals
            local_1=i_x+FNd-1; 
            local_2=i_x+FNd;
            local_3=i_x+2*FNd;
            checkvalue_sigmay=solution_u(local_3);
            if i_x>1
                checkvalue_u=-0.5*(solution_u(local_1)+solution_u(local_2));
            else
                checkvalue_u=-0.5*solution_u(local_2);
            end
            if (state(i_x)==1) && (checkvalue_u<rrr2)
                state(i_x)=2;
                change=1;
            elseif state(i_x)==2 && (checkvalue_u>rrr2)
                state(i_x)=1;
                change=1;
            elseif state(i_x)==2 && (checkvalue_u<rrr1)
                state(i_x)=3;
                change=1;
            elseif state(i_x)==3 && (checkvalue_u>rrr1)
                state(i_x)=2;
                change=1;
            elseif state(i_x)==3 && (checkvalue_u<0)
                state(i_x)=4;
                change=1;  
            elseif state(i_x)==4 && (checkvalue_u>0)
                state(i_x)=3;
                change=1; 
            end
            % end cases
        end
        
        if change==0
            complete=1;
        end
        
        iteration=iteration+1;
    end
    % end active set method loop   %%%%%%%%%%%%%%%%%%%%%%%% 


    U=solution_u(:);  
    figure(3);
    MyShowDisplacement(TRI,U)
    hold on
  
    
    
%end

Ux = U(1:FNd);
Uy = U(FNd+1:2*FNd);

SigmaY=solution_u(2*FNd+1:2*FNd+nx);
% for i=1:intervals
%     if state(i)==1
%         SigmaY(i)=level2;
%     end
% end


TRI.Nodes(TRI.FNodePtrs,:) = TRI.Nodes(TRI.FNodePtrs,:)+[Ux,Uy];  
quiver(TRI.Nodes(2:nx+1,1),TRI.Nodes(2:nx+1,2), zeros(nx,1),-SigmaY ,'r') 

hold off
end