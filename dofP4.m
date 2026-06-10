function [elem2dof,elem2edge,edge,bdDof,freeDof] = dofP4(elem)
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
elem2dof = uint32(zeros(NT,15));
elem2dof(:,1:3) = elem;

%% Two dof on each edge
% edge 1
idx0 = (elem(:,3) > elem(:,2));
elem2dof(idx0,4) = N  + 3*(elem2edge(idx0,1))-2;
elem2dof(idx0,6) = N  + 3*(elem2edge(idx0,1));
elem2dof(~idx0,4)= N  + 3*(elem2edge(~idx0,1));
elem2dof(~idx0,6)= N  + 3*(elem2edge(~idx0,1))-2;
elem2dof(:,5) = N  + 3*(elem2edge(:,1))-1;
% edge 2
idx0 = (elem(:,3) > elem(:,1));
elem2dof(idx0,7) = N  + 3*(elem2edge(idx0,2));
elem2dof(idx0,9) = N  + 3*(elem2edge(idx0,2))-2;
elem2dof(~idx0,7) = N + 3*(elem2edge(~idx0,2))-2;
elem2dof(~idx0,9) = N + 3*(elem2edge(~idx0,2));
elem2dof(:,8) = N + 3*(elem2edge(:,2))-1;

% edge 3
idx0 = (elem(:,2) > elem(:,1));
elem2dof(idx0,10) = N  + 3*(elem2edge(idx0,3))-2;
elem2dof(idx0,12) = N  + 3*(elem2edge(idx0,3));
elem2dof(~idx0,10) = N + 3*(elem2edge(~idx0,3));
elem2dof(~idx0,12) = N  + 3*(elem2edge(~idx0,3))-2;
elem2dof(:,11) = N  + 3*(elem2edge(:,3))-1;



%% Element dof
%elem2dof(:,10:12) = N+2*NE+elem2edge;
elem2dof(:,13:15) = [(N+3*NE+1:N+3*NE+NT)',(N+3*NE+NT+1:N+3*NE+2*NT)',(N+3*NE+2*NT+1:N+3*NE+3*NT)'];
%% Boundary dof
i1(j(3*NT:-1:1)) = 3*NT:-1:1; 
i1 = i1';
bdEdgeIdx = (i1 == i2);
isBdDof = false(N+3*NE+3*NT,1);
isBdDof(edge(bdEdgeIdx,:)) = true;   % boundary node 
idx = find(bdEdgeIdx);
 isBdDof(N+3*idx) = true;      
 isBdDof(N+3*idx-1) = true;
 isBdDof(N+3*idx-2) = true;
% isBdDof(N+3*NE+3*idx) = true;
% isBdDof(N+3*NE+3*idx-1) = true;
% isBdDof(N+3*NE+3*idx-2) = true;
bdDof = find(isBdDof);
freeDof = find(~isBdDof);
end