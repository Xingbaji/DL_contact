function F=MyIntegral(T,G)
 
%�õ����������ɶ����������С�����ε�����������ʼ��M
Nf=length(T.FNodePtrs);
Nt=size(T.Elements,1);


F=zeros(2*Nf,1); %Ҫ�õ���\BigE��\varepsion(u)���������ÿ��С�������Ƿ�Ƭ������

%�õ��ο����������������������ݶȣ�[vs;vt]_{3*2}=[�ݶ�gamma1,�ݶ�gamma2���ݶ�gamma3]
Vs=[-1,1,0];
Vt=[-1,0,1];
qwts=0.5;

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
         F(llr)=F(llr)+ trans.j*qwts * sum( G(i,:).*[Grads1(1,r),0.5*Grads1(2,r),0.5*Grads1(2,r),0]) ;
         F(Nf+llr)=F(Nf+llr)+ trans.j*qwts * sum( G(i,:).*[0,0.5*Grads1(1,r),0.5*Grads1(1,r),Grads1(2,r)]) ;
      end
   end
   
end
 
 
 