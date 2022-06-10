function M=BigE(T,U,coeE,coeK)
% ���M=\mathcal{G} \varepsilon(U)���Ǹ�Nt*4�ľ���

% �õ����������ɶ����������С�����ε�����������ʼ��M
Nf=length(T.FNodePtrs);
Nt=size(T.Elements,1);

coe1=coeE*coeK/(1-coeK^2);
coe2=coeE/(1+coeK);

M=zeros(Nt,4); %������Ҫ�õ���\BigE��\varepsion(u)���������ÿ��С�������Ƿ�Ƭ�����ġ�����

%�õ��ο����������������������ݶȣ�[vs;vt]_{3*2}=[�ݶ�gamma1,�ݶ�gamma2���ݶ�gamma3]
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

