function [coords,ptrs,indices]=MyGetNodes(T,k)

% [coords,ptrs,indices]=getNodes(T,k)
%
%   This function extracts the nodes from triangle k of
%   mesh T.  The output is
%
%       coords: id by 2 array; coordinates of the nodes.
%       ptrs: Pointers into FNodePtrs and CNodePtrs of the
%             nodes (extracted from NodePtrs).
%       indices: Indices of the nodes in Nodes.
%
%   See "help Mesh2" for a description of the mesh T.

%   This routine is part of the MATLAB Fem code that
%   accompanies "Understanding and Implementing the Finite
%   Element Method" by Mark S. Gockenbach (copyright SIAM 2006).

d=T.Degree;
id=(d+1)*(d+2)/2;

% First, get the edges of T_k:

eptrs=T.Elements(k,:);

% Get the indices of the nodes on the edges:

i=T.Edges(abs(eptrs),1:d+1);

% Extract the indices of the vertices and edge nodes:

indices=zeros(id,1);
if eptrs(1)>0
   indices(1:d+1)=i(1,:);
else
   indices(1:d+1)=fliplr(i(1,:));
end
if eptrs(2)>0
   indices(d+1:2*d+1)=i(2,:);
else
   indices(d+1:2*d+1)=fliplr(i(2,:));
end
if eptrs(3)>0
   indices(2*d+1:3*d)=i(3,1:d);
else
   indices(2*d+1:3*d)=fliplr(i(3,2:d+1));
end

% Extract the indices of the interior nodes (if any):

if d>2
   indices(3*d+1:id)=T.IntNodes(k,:);
end

% Extract the coordinates and pointers:

coords=T.Nodes(indices,1:2);
ptrs=T.NodePtrs(indices);