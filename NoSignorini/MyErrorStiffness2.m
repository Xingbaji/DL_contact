function errorE2=MyErrorStiffness2(T,U)

% 得到网格中自由顶点的数量、小三角形的数量，并初始化M
Nf=length(T.FNodePtrs);
Nt=size(T.Elements,1);
MM=sparse(2*Nf,2*Nf);


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

   coeff=trans.j*qwts;

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

               I1=(Grads1(1,r).*Grads1(1,s))*coeff+...
                  (Grads1(2,r).*Grads1(2,s))*coeff/2;
               I2=(Grads1(1,r).*Grads1(1,s))*coeff/2+...
                  (Grads1(2,r).*Grads1(2,s))*coeff;
               I3=(Grads1(2,r).*Grads1(1,s))*coeff/2;
               if r~=s
                  I4=(Grads1(1,r).*Grads1(2,s))*coeff/2;
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
               MM(ii,jj)=MM(ii,jj)+I1;
               MM(ii1,jj1)=MM(ii1,jj1)+I2;
               MM(ii,jj1)=MM(ii,jj1)+I3;
               if r~=s
                  MM(jj,ii1)=MM(jj,ii1)+I4;
               end
            end
         end
      end
   end
end

% Fill in the lower triangle of K, using symmetry:
MM=MM+triu(MM,1)';


Umatrix=U'*U;
MM=MM.*Umatrix;

errorE2=sqrt(sum(MM(:)));

