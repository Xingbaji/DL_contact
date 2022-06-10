function trans=MyTransToRefTri(c)

% trans=TransToRefTri(c)
%
%   This function takes as input a 3 by 2 array c defining
%   the vertices of a triangle T, and computes the
%   transformation from the reference triangle with vertices
%   (0,0), (1,0), (0,1) to T.
%
%   The output trans is a struct with the following fields:
%
%      trans.J - the Jacobian matrix of the transformation
%      trans.j - the absolute value of det(trans.J)
%      trans.z1- the first vertex of the triangle

%   This routine is part of the MATLAB Fem code that
%   accompanies "Understanding and Implementing the Finite
%   Element Method" by Mark S. Gockenbach (copyright SIAM 2006).

% Compute the Jacobian matrix and its determinant

trans.J=[
c(2,1)-c(1,1) c(3,1)-c(1,1)
c(2,2)-c(1,2) c(3,2)-c(1,2)];

trans.j=abs(det(trans.J));

% Record the first vertex:

trans.z1=c(1,:)';
