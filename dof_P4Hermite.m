function [elem2dof,elem2edge,edge,bdDof,freeDof] = dof_P4Hermite(elem,node)
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
elem2dof = [elem,elem + N,elem+2*N,3*N+elem2edge,(3*N+NE+1:3*N+NE+NT)'....
    ,(3*N+NE+NT+1:3*N+NE+2*NT)',(3*N+NE+2*NT+1:3*N+NE+3*NT)'];


%% Element dof
%elem2dof(:,10:12) = N+2*NE+elem2edge;
%% Boundary dof
i1(j(3*NT:-1:1)) = 3*NT:-1:1; 
i1 = i1';
bdEdgeIdx = (i1 == i2);
isBdDof = false(3*N+NE+3*NT,1);
isBdDof(edge(bdEdgeIdx,:)) = true;   % boundary node 
bdNode = double(unique(edge(bdEdgeIdx,:)));
bx1 = find(node(:,1) == 0);
bx2 = find(node(:,1) == 1);
by1 = find(node(:,2) == 0);
by2 = find(node(:,2) == 1);
isBdDof(N + [by1;by2]) = true;
isBdDof(2*N + [bx1;bx2]) = true;
idx = find(bdEdgeIdx);
isBdDof(3*N+idx) = true;  
%isBdDof(bdNode+3*N) = ture;
%isBdDof(bdNode+4*N) = true; % Dxy
%bx1 = find(node(:,1) == 0);
%bx2 = find(node(:,1) == 1);
%by1 = find(node(:,2) == 0);
%by2 = find(node(:,2) == 1);
%isBdDof(3*N + [bx1;bx2]) = true;
%isBdDof(5*N + [by1;by2]) = true;
% zuo = find(node(:,1) == -1);
% you = find(node(:,1) ==  1);
% shang = find(node(:,2)== 1);
% xia = find(node(:,2)== -1);
% isBdDof(xia+N) = true;
% isBdDof(shang+N) = true;
% isBdDof(zuo+2*N) = true;
% isBdDof(you+2*N) = true;

bdDof = find(isBdDof);
freeDof = find(~isBdDof);
end