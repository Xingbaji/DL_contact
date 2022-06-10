function MyShowDisplacement(T,U,g)

% ShowDisplacement(T,U,g)
%
%   This function displays a mesh and the corresponding
%   displaced mesh, where the nodal values U (free nodes)
%   and g (constrained nodes) define the displacement.
%
%   See "help Mesh" for a description of the mesh T.
%
%   U is a (2*Nf) vector, where Nf is the number of free
%   nodes.  The first Nf components of U contain the
%   horizontal components of the displacements at the
%   free nodes, and the second Nf components contain
%   the vertical components.
%
%   g is a (2*Nc) vector, where Nc is the number of
%   constrained nodes, and is formed similarly.

%   This routine is part of the MATLAB Fem code that
%   accompanies "Understanding and Implementing the Finite
%   Element Method" by Mark S. Gockenbach (copyright SIAM 2006).

Nf=length(T.FNodePtrs);
Nc=length(T.CNodePtrs);
if nargin<3||isempty(g)
   g=zeros(2*Nc,1);
end

MyShowMesh(T,[0,0,0,0],'k')
hold on
T.Nodes(T.FNodePtrs,:)=T.Nodes(T.FNodePtrs,:)+[U(1:Nf),U(Nf+1:2*Nf)];
if Nc>0
   T.Nodes(T.CNodePtrs,:)=T.Nodes(T.CNodePtrs,:)+[g(1:Nc),g(Nc+1:2*Nc)];
end
MyShowMesh(T,[0,0,0,0],'r')
hold off
