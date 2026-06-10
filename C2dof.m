function [elem2dof,elem2edge,edge,bdDof,freeDof] = C2dof(elem,node)
%% DOFP3 dof structure for P3 element.
%
%  [elem2dof,edge,bdDof] = DOFP3(elem) constructs the dof structure
%  for the quadratic element based on a triangle. elem2dof(t,i) is the
%  global index of the i-th dof of the t-th element.
%
%  The global indices of the dof is organized  according to the order of
%  nodes, edges and elements, namely, first give index number to the dofs
%  on nodes, then the dofs on edges, last the dofs on elements.
%
%  See also dofP2, dof3P3.
%  
%  Created by Jie Zhou.  
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

N = max(max(elem)); NT = size(elem,1);  

%% Data structure
totalEdge = uint32(sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2));
[edge, i2, j] = myunique(totalEdge);
NE = size(edge,1);
elem2edge = reshape(j,NT,3);

%% Nodal dof
%elem2dof = uint32(zeros(NT,21));
elem2dof = [elem,elem + N,elem+2*N,elem+3*N,elem+4*N,elem+5*N,(6*N+1:6*N+NT)',(6*N+NT+1:6*N+2*NT)',(6*N+2*NT+1:6*N+3*NT)'];


%% Element dof
%elem2dof(:,10:12) = N+2*NE+elem2edge;
%% Boundary dof
i1(j(3*NT:-1:1)) = 3*NT:-1:1; 
i1 = i1';
bdEdgeIdx = (i1 == i2);
isBdDof = false(6*N+3*NT,1);
isBdDof(edge(bdEdgeIdx,:)) = true;   % boundary node 
bdNode = unique(edge(bdEdgeIdx,:));
%isBdDof(bdNode+N) = true;
%isBdDof(bdNode+2*N) = true;
%isBdDof(bdNode+3*N) = ture;
%isBdDof(bdNode+4*N) = true; % Dxy
by1 = find(node(:,1) == -1);
by2 = find(node(:,1) == 1);
%by3 = find(node(:,1) == 0 & node(:,2)<= 0);
bx1 = find(node(:,2) == -1);
bx2 = find(node(:,2) == 1);
%bx3 = find(node(:,2) == 0 & node(:,1)>=0);
isBdDof(N + [bx1;bx2]) = true;
isBdDof(2*N + [by1;by2]) = true;
isBdDof(3*N + [bx1;bx2]) = true;
isBdDof(5*N + [by1;by2]) = true;
%idx = find(bdEdgeIdx);
% isBdDof(6*N+idx) = true; 
bdDof = find(isBdDof);
freeDof = find(~isBdDof);
end