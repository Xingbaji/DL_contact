function M=BigE(T,U,coeE,coeK)
% 输出M=\mathcal{G} \varepsilon(U)，是个Nt*4的矩阵

% 得到网格中自由顶点的数量、小三角形的数量，并初始化M
Nf=length(T.FNodePtrs);
Nt=size(T.Elements,1);

coe1=coeE*coeK/(1-coeK^2);
coe2=coeE/(1+coeK);

M=zeros(Nt,4); %《《《要得到的\BigE（\varepsion(u)）在网格的每个小三角形是分片常数的》》》

%得到参考三角形上三个基函数的梯度：[vs;vt]_{3*2}=[梯度gamma1,梯度gamma2，梯度gamma3]
Vs=[-1,1,0];
Vt=[-1,0,1];

% Add the contributions from each element
for i=1:Nt
    
   % Get the coordinates and pointers of the nodes:
   [coords,ll]=MyGetNodes(T,i);

   % Transform the reference triangle to T. 
   trans=MyTransToRefTri(coords);

   Grads1=trans.J'\[Vs;Vt];

   % Compute the contributions to K by looping over all possible
   % combinations of (global) indices related to this triangle.
   for r=1:3
      llr=ll(r);
      if llr>0
         M(i,1)=M(i,1)+ U(llr)*(coe1+coe2)*Grads1(1,r) +U(Nf+llr)*coe1*Grads1(2,r) ;
         M(i,2)=M(i,2)+ U(llr)*coe2/2*Grads1(2,r) +U(Nf+llr)*coe2/2*Grads1(1,r) ;
         M(i,3)=M(i,3)+ U(llr)*coe2/2*Grads1(2,r) +U(Nf+llr)*coe2/2*Grads1(1,r) ;
         M(i,4)=M(i,4)+ U(llr)*coe1*Grads1(1,r) +U(Nf+llr)*(coe1+coe2)*Grads1(2,r) ;
      end
   end
end

