function K=MyStiffness2(T,coeE,coeK)
%经2*2网格笔算，证明此函数正确

%  参考 "Understanding and Implementing the Finite Element Method" by Mark S.
%  Gockenbach (copyright SIAM 2006)程序包version3中的“function
%  K=StiffnessIso(T,mu,lam)”.

Coe1=coeE/(1+coeK);
Coe2=coeE*coeK/(1-coeK^2);

% 得到网格中自由顶点的数量、小三角形的数量，并初始化K
Nf=length(T.FNodePtrs);
Nt=size(T.Elements,1);
K=sparse(2*Nf,2*Nf);

%参考三角形的面积是1/2
qwts=1/2;   

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


   % Compute the values of mu and lam at the quadrature
   % nodes and multiply them by the quadrature weights
   % and Jacobians for convenience:
   muvals=Coe1*qwts;
   lamvals=Coe2*qwts;

   lamvals=trans.j*lamvals;
   muvals=trans.j*muvals;
   mulamvals=muvals+lamvals;

   % Compute the contributions to K by looping over all possible
   % combinations of (global) indices related to this triangle.

   for r=1:3
      llr=ll(r);
      if llr>0
         for s=r:3
            lls=ll(s);
            if lls>0

               % These basis functions on this element contribute
               % to four entries in the upper triangle of K
               % (three if ll(r)==ll(s)).

               I1=(Grads1(1,r).*Grads1(1,s))*mulamvals+...
                  (Grads1(2,r).*Grads1(2,s))*muvals/2;
               I2=(Grads1(1,r).*Grads1(1,s))*muvals/2+...
                  (Grads1(2,r).*Grads1(2,s))*mulamvals;
               I3=(Grads1(2,r).*Grads1(1,s))*muvals/2+...
                  (Grads1(1,r).*Grads1(2,s))*lamvals;
               if r~=s
                  I4=(Grads1(1,r).*Grads1(2,s))*muvals/2+...
                     (Grads1(2,r).*Grads1(1,s))*lamvals;
               end

               if llr>lls
                  ii=lls;
                  jj=llr;
                  tmp=I3;
                  I3=I4;
                  I4=tmp;
               else
                  ii=llr;
                  jj=lls;
               end
               ii1=ii+Nf;
               jj1=jj+Nf;
               K(ii,jj)=K(ii,jj)+I1;
               K(ii1,jj1)=K(ii1,jj1)+I2;
               K(ii,jj1)=K(ii,jj1)+I3;
               if r~=s
                  K(jj,ii1)=K(jj,ii1)+I4;
               end
            end
         end
      end
   end
end

% Fill in the lower triangle of K, using symmetry:
K=K+triu(K,1)';