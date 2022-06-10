function MyShowMesh(T,flags,c)

% ShowMesh(T,flags,c)
%
%   This function displays a triangular mesh T.  For a
%   description of the data structure T, see "help Mesh".
%
%   The optional argument flags is an array of flags
%   with the following effects:
%
%    flags(1)==1: the nodes are labeled by their indices
%              2: the nodes are labeled by their NodePtrs
%    flags(2): the edges are labeled by their indices
%    flags(3): the triangles are labeled by their indices
%    flags(4)==1: the nodes are indicated by '.'
%              2: the nodes are indicated by 'o' for free and
%                 a star for constrained
%
%   (flags(2) and flags(3) are simply zero (off) or
%   nonzero (on).)
%
%   The optional input c is the color ('b', 'r', etc) of
%   the edges.

%   This routine is part of the MATLAB Fem code that
%   accompanies "Understanding and Implementing the Finite
%   Element Method" by Mark S. Gockenbach (copyright SIAM 2006).

% Assign the optional arguments, if necessary:

if nargin<3
   c='k';
end

if nargin<2
   flags=zeros(1,4);
elseif length(flags)<4
   flags=[flags(:)',zeros(1,4-length(flags))];
end

if flags(1)>0 && flags(4)==0
   flags(4)=1;
end

d=T.Degree;

% Display the true triangles:

X=T.Nodes(:,1);
Y=T.Nodes(:,2);
[tris,tris1]=MygetTriNodeIndices3(T);
tris=tris(:,1:d:2*d+1);
trimesh(tris,X,Y,zeros(size(X)),'facecolor','none','edgecolor',c);
hold on

% Display the triangles with curved edges:

Nt1=size(tris1,1);
if Nt1>0

   TR=RefTri(d);
   inodes=getNodes(TR,1);
   s=linspace(0,1,21)';
   enodes=[s,1-s];
   Vals=EvalNodalBasisFcns(inodes,enodes);
   for i=1:Nt1

      nodes=T.Nodes(tris1(i,:),:);
      nodes1=Vals*nodes;
      x=[nodes(1,1);nodes1(:,1);nodes(1,1)];
      y=[nodes(1,2);nodes1(:,2);nodes(1,2)];
      plot(x,y,c,'LineWidth',1)

   end

end

if flags(4)==1
   hold on
   plot3(X,Y,zeros(size(X)),'.','MarkerSize',10)
end
view(2)
axis('equal','off')
if any(flags)

   tris=MygetTriNodeIndices(T);

   % Extract the number of triangles, edges, and nodes

   Nt=size(T.Elements,1);
   Ne=size(T.Edges,1);
   Nv=length(T.Nodes);

   % Label the triangles, if desired.

   if flags(3)

      for i=1:Nt
         x=mean(X(tris(i,:)));
         y=mean(Y(tris(i,:)));
         text(x,y,int2str(i),'FontSize',12,'Color','black')
      end

   end

   % Display the nodes and label them, if desired.

   del=(max(abs(T.Nodes(:,1)))-min(abs(T.Nodes(:,1))))/50;

   if flags(1) || flags(4)==2

      for i=1:Nv

         if T.NodePtrs(i)>0 && flags(4)==2
            plot(T.Nodes(i,1),T.Nodes(i,2),['o','b'])
         elseif T.NodePtrs(i)<0 && flags(4)==2
            plot(T.Nodes(i,1),T.Nodes(i,2),['h','b'])
         end

         if flags(1)==1
            text(T.Nodes(i,1)+del,T.Nodes(i,2),int2str(i),...
                 'FontSize',12,'Color','blue')
         end

         if flags(1)==2
            text(T.Nodes(i,1)+del,T.Nodes(i,2),int2str(T.NodePtrs(i)),...
                 'FontSize',12,'Color','blue')
         end

      end

   end

   % Label the edges, if desired.

   if flags(2)

      for i=1:Ne
         c=0.5*sum(T.Nodes(T.Edges(i,[1,d+1]),:));
         text(c(1),c(2),int2str(i),'FontSize',12,'Color','red')
      end

   end

end

hold off

return

function [ElList,ElList1]=MygetTriNodeIndices3(T)

% [ElList,ElList1]=getTriNodeIndices3(T)
%
%   This function creates an Nt1 by id=(d+1)(d+2)/2 array, and
%   an Nt2 by id array.  Each row corresponds to one triangle
%   in the mesh T and contains the indices of the nodes in
%   T.Nodes (d is the degree of the triangles).
%
%   The array ElList contains the true triangles, and
%   the array ElList1 contains those "triangles" having a
%   curved edge.
%
%   For a description of the data structure T, see "help Mesh".

d=T.Degree;
id=(d+1)*(d+2)/2;

Nt=size(T.Elements,1);
ElList=zeros(Nt,id);
ElList1=zeros(Nt,id);
n1=0;
n2=0;

for k=1:Nt
   [t1,t2,indices]=MyGetNodes(T,k);
   if any(T.EdgeCFlags(abs(T.Elements(k,:)))>0)
      n2=n2+1;
      ElList1(n2,:)=indices';
   else
      n1=n1+1;
      ElList(n1,:)=indices';
   end
end
ElList=ElList(1:n1,:);
ElList1=ElList1(1:n2,:);
