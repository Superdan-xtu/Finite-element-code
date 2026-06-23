 
%% Introduction
% This script runs the verification experiment for the Hermite-type
% finite element pair for the linear stress gradient elasticity problem.
%
% The script can be run directly in MATLAB. It computes numerical
% solutions on a sequence of meshes and compares them with a precomputed
% reference solution obtained on a fine mesh.
%
% Current setting:
%   - finite element: P4 Hermite-type stress element
%   - right-hand side: La_Ipf_LSG.m
%   - reference solution: solutiota10e4new.mat
%   - small parameter: iota = 1.0e-04
%
% Important:
% The value of iota must be consistent with the loaded reference data.
% For example,
%
%   iota = 1          -> solutiota1new.mat
%   iota = 1.0e-01    -> solutiota10e1new.mat
%   iota = 1.0e-02    -> solutiota10e2new.mat
%   iota = 1.0e-04    -> solutiota10e4new.mat
%
% To test another value of iota, change both the line defining iota and
% the corresponding reference solution file loaded below.
%
% Output:
% The script prints the mesh size, the L2 error of the displacement, and
% the H1-type error of the stress in the MATLAB command window.

for mesh_s=1:1
%%  Mesh Setting
% hsize = 1/(2*mesh_s);    
% [node1,elem1] = squaremesh([0,1,0,1],hsize); 

La = 1.0;
mu = 0.3;
iota = 1.0;
%[RDlambda,Rarea,RelemSign] = gradbasis(Rnode,Relem);
%[~,~,Relem2edgeSign,RedgeSign,Redge2elem,Rbdedge] = dofedgemy(Relem);
%[ elemt , elemn ] = elemtau_n( node , Relem , Relem2edgeSign ,Relem2edge ,Redge );
%[Relem2dof1,Relem2edge,Redge,RbdDof,RfreeDof] = dofP4(Relem); 
%[Relem2dof1,Relem2edge,Redge,RbdDof,RfreeDof] = C2dof(Relem,Rnode); 
%Relem2dof = double(Relem2dof1);


base_num = 36;
% inter_node2 = 0.4*node1(elem1(:,1),:)+0.2*node1(elem1(:,2),:)+0.4*node1(elem1(:,3),:);
% inter_node3 = 0.4*node1(elem1(:,1),:)+0.4*node1(elem1(:,2),:)+0.2*node1(elem1(:,3),:);
% inter_node1 = 0.2*node1(elem1(:,1),:)+0.4*node1(elem1(:,2),:)+0.4*node1(elem1(:,3),:);
% N = size(node1,1);
% NT = size(elem1,1);
% node = [node1;inter_node1;inter_node2;inter_node3];
% elem_inter1 = [(N+1:N+NT)',(N+NT+1:N+2*NT)',(N+2*NT+1:N+3*NT)'];
% elem_e1 = [elem1(:,2),elem1(:,3),(N+1:N+NT)' ];
% elem_e2 = [(N+NT+1:N+2*NT)',elem1(:,3),elem1(:,1)];
% elem_e3 = [(N+2*NT+1:N+3*NT)',elem1(:,1),elem1(:,2)];
% elem_v1 = [ (N + NT + 1: N+ 2*NT)',elem1(:,1), (N+2*NT+1:N+3*NT)' ];
% elem_v2 = [ (N + 1:N+NT)', (N+2*NT+1:N+3*NT)', elem1(:,2) ];
% elem_v3 = [ (N + 1:N+NT)', elem1(:,3), (N+NT+1:N+2*NT)'  ];
% 
% elem = [elem_inter1;elem_v1;elem_v2;elem_v3;elem_e1;elem_e2;elem_e3];

 % inter_node1 = 1/3*node1(elem1(:,1),:) + 1/3*node1(elem1(:,2),:) + 1/3*node1(elem1(:,3),:);
 %e_node =  1/2*node1(edge1(:,2),:) + 1/2*node1(edge1(:,1),:);
 %e2_node =  1/2*node1(elem1(:,1),:) + 1/2*node1(elem1(:,3),:);
 %e3_node =  1/2*node1(elem1(:,2),:) + 1/2*node1(elem1(:,1),:);
% node = [node1;inter_node1];
% N = size(node1,1);
% NT = size(elem1,1);
%elem_v1 = [elem1(:,2),(N+1:N+NT)',elem1(:,3)];
%elem_v2 = [elem1(:,1),(N+1:N+NT)',elem1(:,3)]; 
%elem_v3 = [elem1(:,1),(N+1:N+NT)',elem1(:,2)];

%elem = [elem_v1;elem_v2;elem_v3]; 
delta = 1.0e-4;%0.0e-01;
x_node = 1/2 + delta;
y_node = 1/2 - delta;
node = [0,0;1,0;1,1;0,1;x_node,y_node];
elem = [1,2,5;2,3,5;5,3,4;4,1,5];
% 
% elem = double(elem);
%node = [0,0;1,0;0,1];
%elem = [1,2,3];

%Mesh_G = size(Relem,1)/size(elem,1);

[Dlambda,area,elemSign] = gradbasis(node,elem);
[elem2edge,edge,elem2edgeSign,edgeSign,edge2elem,bdedge] = dofedgemy(elem);
[ elemt , elemn ] = elemtau_n( node , elem , elem2edgeSign ,elem2edge ,edge );
elem = double(elem);
[elem2dof1,elem2edge,edge,bdDof,freeDof] = dof_P7Hermite(elem,node); 
Ndof = double(max(elem2dof1(:)));
elem2dof = [elem2dof1];
Stiff = sparse(Ndof,Ndof);
M =  Stiff;
% ç»Ÿä¸€æ³•å‘
elemnor = elemn;
elemn(:,:,1)=elemn(:,:,1).*elem2edgeSign(:,1);
elemn(:,:,2)=elemn(:,:,2).*elem2edgeSign(:,2);
elemn(:,:,3)=elemn(:,:,3).*elem2edgeSign(:,3);

order1 = 9;
[lambda,weight1] = quadpts1(order1);
nQuad = size(lambda,1);
order2 = 20;
[lambda2,weight2] = quadptsmy(order2);
Quad2= size(lambda2,1);

test_q = 0;
 for p = 1:Quad2
 test_q = test_q + weight2(p)*lambda2(p,1)^15;
 end
DG_base_num = 28;
DG_phi = zeros(Quad2,DG_base_num);

DG_phi(:,1)  = lambda2(:,1).^6;
DG_phi(:,2)  = lambda2(:,2).^6;
DG_phi(:,3)  = lambda2(:,3).^6;

DG_phi(:,4)  = lambda2(:,1).^5.*lambda2(:,2);
DG_phi(:,5)  = lambda2(:,1).^5.*lambda2(:,3);
DG_phi(:,6)  = lambda2(:,2).^5.*lambda2(:,1);
DG_phi(:,7)  = lambda2(:,2).^5.*lambda2(:,3);
DG_phi(:,8)  = lambda2(:,3).^5.*lambda2(:,1);
DG_phi(:,9)  = lambda2(:,3).^5.*lambda2(:,2);

DG_phi(:,10) = lambda2(:,1).^4.*lambda2(:,2).^2;
DG_phi(:,11) = lambda2(:,1).^4.*lambda2(:,3).^2;
DG_phi(:,12) = lambda2(:,2).^4.*lambda2(:,1).^2;
DG_phi(:,13) = lambda2(:,2).^4.*lambda2(:,3).^2;
DG_phi(:,14) = lambda2(:,3).^4.*lambda2(:,1).^2;
DG_phi(:,15) = lambda2(:,3).^4.*lambda2(:,2).^2;

DG_phi(:,16) = lambda2(:,1).^4.*lambda2(:,2).*lambda2(:,3);
DG_phi(:,17) = lambda2(:,2).^4.*lambda2(:,1).*lambda2(:,3);
DG_phi(:,18) = lambda2(:,3).^4.*lambda2(:,1).*lambda2(:,2);

DG_phi(:,19) = lambda2(:,1).^3.*lambda2(:,2).^3;
DG_phi(:,20) = lambda2(:,1).^3.*lambda2(:,3).^3;
DG_phi(:,21) = lambda2(:,2).^3.*lambda2(:,3).^3;

DG_phi(:,22) = lambda2(:,1).^3.*lambda2(:,2).^2.*lambda2(:,3);
DG_phi(:,23) = lambda2(:,1).^3.*lambda2(:,2).*lambda2(:,3).^2;
DG_phi(:,24) = lambda2(:,2).^3.*lambda2(:,1).^2.*lambda2(:,3);
DG_phi(:,25) = lambda2(:,2).^3.*lambda2(:,1).*lambda2(:,3).^2;
DG_phi(:,26) = lambda2(:,3).^3.*lambda2(:,1).^2.*lambda2(:,2);
DG_phi(:,27) = lambda2(:,3).^3.*lambda2(:,1).*lambda2(:,2).^2;

DG_phi(:,28) = lambda2(:,1).^2.*lambda2(:,2).^2.*lambda2(:,3).^2;
N = size(node,1);
NT = size(elem,1);
DG_dof = [elem, ...
    (N+1:N+NT)', ...
    (N+NT+1:N+2*NT)', ...
    (N+2*NT+1:N+3*NT)', ...
    (N+3*NT+1:N+4*NT)', ...
    (N+4*NT+1:N+5*NT)', ...
    (N+5*NT+1:N+6*NT)', ...
    (N+6*NT+1:N+7*NT)', ...
    (N+7*NT+1:N+8*NT)', ...
    (N+8*NT+1:N+9*NT)', ...
    (N+9*NT+1:N+10*NT)', ...
    (N+10*NT+1:N+11*NT)', ...
    (N+11*NT+1:N+12*NT)', ...
    (N+12*NT+1:N+13*NT)', ...
    (N+13*NT+1:N+14*NT)', ...
    (N+14*NT+1:N+15*NT)', ...
    (N+15*NT+1:N+16*NT)', ...
    (N+16*NT+1:N+17*NT)', ...
    (N+17*NT+1:N+18*NT)', ...
    (N+18*NT+1:N+19*NT)', ...
    (N+19*NT+1:N+20*NT)', ...
    (N+20*NT+1:N+21*NT)', ...
    (N+21*NT+1:N+22*NT)', ...
    (N+22*NT+1:N+23*NT)', ...
    (N+23*NT+1:N+24*NT)', ...
    (N+24*NT+1:N+25*NT)'];
NDG = max(max(DG_dof));

BStiff = sparse(3*Ndof,2*NDG);
M_DG = sparse(NDG,NDG);
DG_base_num = 28;

DGbase = zeros(DG_base_num,DG_base_num);
DGbase(1,1) = 1;
DGbase(2,2) = 1;
DGbase(3,3) = 1;

DGm = [...
    6,0,0;
    0,6,0;
    0,0,6;
    5,1,0;
    5,0,1;
    1,5,0;
    0,5,1;
    1,0,5;
    0,1,5;
    4,2,0;
    4,0,2;
    2,4,0;
    0,4,2;
    2,0,4;
    0,2,4;
    4,1,1;
    1,4,1;
    1,1,4;
    3,3,0;
    3,0,3;
    0,3,3;
    3,2,1;
    3,1,2;
    2,3,1;
    1,3,2;
    2,1,3;
    1,2,3;
    2,2,2];

for i = 1:DG_base_num
    for j = 4:DG_base_num
        DGbase(j,i) = com_quad((DGm(j,:)+DGm(i,:))');
    end
end

DGcof = inv(DGbase);
%DGcof = DGcof'; %% ???????????????
phi = zeros(Quad2,base_num);

phi(:,1)  = lambda2(:,1).^7;
phi(:,2)  = lambda2(:,2).^7;
phi(:,3)  = lambda2(:,3).^7;

phi(:,4)  = lambda2(:,1).^6.*lambda2(:,2);
phi(:,5)  = lambda2(:,1).^6.*lambda2(:,3);
phi(:,6)  = lambda2(:,2).^6.*lambda2(:,1);
phi(:,7)  = lambda2(:,2).^6.*lambda2(:,3);
phi(:,8)  = lambda2(:,3).^6.*lambda2(:,1);
phi(:,9)  = lambda2(:,3).^6.*lambda2(:,2);

phi(:,10) = lambda2(:,1).^5.*lambda2(:,2).^2;
phi(:,11) = lambda2(:,1).^5.*lambda2(:,3).^2;
phi(:,12) = lambda2(:,2).^5.*lambda2(:,1).^2;
phi(:,13) = lambda2(:,2).^5.*lambda2(:,3).^2;
phi(:,14) = lambda2(:,3).^5.*lambda2(:,1).^2;
phi(:,15) = lambda2(:,3).^5.*lambda2(:,2).^2;

phi(:,16) = lambda2(:,1).^5.*lambda2(:,2).*lambda2(:,3);
phi(:,17) = lambda2(:,2).^5.*lambda2(:,1).*lambda2(:,3);
phi(:,18) = lambda2(:,3).^5.*lambda2(:,1).*lambda2(:,2);

phi(:,19) = lambda2(:,1).^4.*lambda2(:,2).^3;
phi(:,20) = lambda2(:,1).^4.*lambda2(:,3).^3;
phi(:,21) = lambda2(:,2).^4.*lambda2(:,1).^3;
phi(:,22) = lambda2(:,2).^4.*lambda2(:,3).^3;
phi(:,23) = lambda2(:,3).^4.*lambda2(:,1).^3;
phi(:,24) = lambda2(:,3).^4.*lambda2(:,2).^3;

phi(:,25) = lambda2(:,1).^4.*lambda2(:,2).^2.*lambda2(:,3);
phi(:,26) = lambda2(:,1).^4.*lambda2(:,2).*lambda2(:,3).^2;
phi(:,27) = lambda2(:,2).^4.*lambda2(:,1).^2.*lambda2(:,3);
phi(:,28) = lambda2(:,2).^4.*lambda2(:,1).*lambda2(:,3).^2;
phi(:,29) = lambda2(:,3).^4.*lambda2(:,1).^2.*lambda2(:,2);
phi(:,30) = lambda2(:,3).^4.*lambda2(:,1).*lambda2(:,2).^2;

phi(:,31) = lambda2(:,1).^3.*lambda2(:,2).^3.*lambda2(:,3);
phi(:,32) = lambda2(:,1).^3.*lambda2(:,2).*lambda2(:,3).^3;
phi(:,33) = lambda2(:,1).*lambda2(:,2).^3.*lambda2(:,3).^3;

phi(:,34) = lambda2(:,1).^3.*lambda2(:,2).^2.*lambda2(:,3).^2;
phi(:,35) = lambda2(:,2).^3.*lambda2(:,1).^2.*lambda2(:,3).^2;
phi(:,36) = lambda2(:,3).^3.*lambda2(:,1).^2.*lambda2(:,2).^2;

% P7 basis exponents, consistent with the 36 basis functions phi(:,1:36)
m7 = [...
    7,0,0;
    0,7,0;
    0,0,7;
    6,1,0;
    6,0,1;
    1,6,0;
    0,6,1;
    1,0,6;
    0,1,6;
    5,2,0;
    5,0,2;
    2,5,0;
    0,5,2;
    2,0,5;
    0,2,5;
    5,1,1;
    1,5,1;
    1,1,5;
    4,3,0;
    4,0,3;
    3,4,0;
    0,4,3;
    3,0,4;
    0,3,4;
    4,2,1;
    4,1,2;
    2,4,1;
    1,4,2;
    2,1,4;
    1,2,4;
    3,3,1;
    3,1,3;
    1,3,3;
    3,2,2;
    2,3,2;
    2,2,3];

base_num = 36;

%% interior moments: P4
% test basis exponents for P4
m4 = [...
    4,0,0;
    0,4,0;
    0,0,4;
    3,1,0;
    3,0,1;
    1,3,0;
    0,3,1;
    1,0,3;
    0,1,3;
    2,2,0;
    2,0,2;
    0,2,2;
    2,1,1;
    1,2,1;
    1,1,2];

BA = zeros(15,base_num);
for j = 1:15
    for i = 1:base_num
        BA(j,i) = com_quad((m7(i,:) + m4(j,:))');
    end
end

%% edge moments: P3 on each edge
% 1D P3 basis on an edge, written by barycentric powers on the two nonzero lambdas
mE3 = [...
    3,0;
    2,1;
    1,2;
    0,3];

BAE = zeros(12,base_num);

for i = 1:base_num
    % edge 1: lambda1 = 0, use (lambda2, lambda3)
    for j = 1:4
        BAE(j,i) = 0^m7(i,1) * com_quad((m7(i,[2,3]) + mE3(j,:))');
    end

    % edge 2: lambda2 = 0, use (lambda1, lambda3)
    for j = 1:4
        BAE(4+j,i) = 0^m7(i,2) * com_quad((m7(i,[3,1]) + mE3(j,:))');
    end

    % edge 3: lambda3 = 0, use (lambda1, lambda2)
    for j = 1:4
        BAE(8+j,i) = 0^m7(i,3) * com_quad((m7(i,[1,2]) + mE3(j,:))');
    end
end

% for i = 1:base_num
%     % edge 1: lambda1 = 0, use (lambda2, lambda3)
%     for j = 1:4
%         BAE(3*j-2,i) = 0^m7(i,1) * com_quad((m7(i,[2,3]) + mE3(j,:))');
%         BAE(3*j-1,i) = 0^m7(i,2) * com_quad((m7(i,[3,1]) + mE3(j,:))');  
%         BAE(3*j,i)  =  0^m7(i,3) * com_quad((m7(i,[1,2]) + mE3(j,:))');
%     end
% end

COF = zeros(base_num*NT,base_num);



for i = 1:NT
nor = zeros(3,2);
nor(1,:)=elemn(i,:,1);
nor(2,:)=elemn(i,:,2);
nor(3,:)=elemn(i,:,3);
Dl = zeros(3,2);
Dl(1,:)=Dlambda(i,:,1);
Dl(2,:)=Dlambda(i,:,2);
Dl(3,:)=Dlambda(i,:,3);
nor = nor';
Dl=Dl';

[ A] = combase_P7Hermite( Dl,base_num);
A(10:21,:) = BAE; 
A(22:36,:) = BA;

cof = inv(A);
COF((i-1)*base_num+1:i*base_num,:)=cof'; %æ¯ä¸€è¡Œä»£è¡¨ä¸€ä¸ªå•å…ƒçš„ä¸?ä¸ªåŸºçš„ç³»æ•?

end
clear cof BA A BAE
%% Cof on Refine mesh







AA = zeros(NT,(base_num*(base_num+1)/2));
MM = AA;


elem2dof = double(elem2dof);
for p = 1:Quad2
Dx = zeros(NT,base_num);
Dy = zeros(NT,base_num);

l1 = lambda2(p,1);
l2 = lambda2(p,2);
l3 = lambda2(p,3);

for i = 1:base_num
    a = m7(i,1);
    b = m7(i,2);
    c = m7(i,3);

    if a > 0
        Dx(:,i) = Dx(:,i) + a * l1.^(a-1) .* l2.^b .* l3.^c * Dlambda(:,1,1);
        Dy(:,i) = Dy(:,i) + a * l1.^(a-1) .* l2.^b .* l3.^c * Dlambda(:,2,1);
    end
    if b > 0
        Dx(:,i) = Dx(:,i) + b * l1.^a .* l2.^(b-1) .* l3.^c * Dlambda(:,1,2);
        Dy(:,i) = Dy(:,i) + b * l1.^a .* l2.^(b-1) .* l3.^c * Dlambda(:,2,2);
    end
    if c > 0
        Dx(:,i) = Dx(:,i) + c * l1.^a .* l2.^b .* l3.^(c-1) * Dlambda(:,1,3);
        Dy(:,i) = Dy(:,i) + c * l1.^a .* l2.^b .* l3.^(c-1) * Dlambda(:,2,3);
    end
end

    
for i = 1:base_num
    for j = 1:i
    AA(:,j+(i-1)*i/2) = AA(:,j+(i-1)*i/2) + weight2(p)*(dot(COF(i:base_num:end,:),Dx,2).*dot(COF(j:base_num:end,:),Dx,2)+dot(COF(i:base_num:end,:),Dy,2).*dot(COF(j:base_num:end,:),Dy,2));
    MM(:,j+(i-1)*i/2) = MM(:,j+(i-1)*i/2) + weight2(p)*(COF(i:base_num:end,:)*(phi(p,:)').*COF(j:base_num:end,:)*(phi(p,:)'));
   
  %  Stiff = Stiff + sparse(double(elem2dof(:,i)),double(elem2dof(:,j)),Aij,Ndof,Ndof);    
        
    end
end    
    
    
end
for k = 1:(base_num*(base_num+1)/2)
AA(:,k)= AA(:,k).*area;
MM(:,k) = MM(:,k).*area;
end
for i = 1:base_num
    for j = 1:i
        Aij = AA(:,j+(i-1)*i/2);
        Mij = MM(:,j+(i-1)*i/2);
  if j==i
            Stiff = Stiff + sparse(elem2dof(:,i),elem2dof(:,j),Aij,Ndof,Ndof);
            M = M + sparse(elem2dof(:,i),elem2dof(:,j),Mij,Ndof,Ndof);
  else
            Stiff = Stiff + sparse([elem2dof(:,i);elem2dof(:,j)],[elem2dof(:,j);elem2dof(:,i)],...
                           [Aij; Aij],Ndof,Ndof);       
              M =     M +  sparse([elem2dof(:,i);elem2dof(:,j)],[elem2dof(:,j);elem2dof(:,i)],...
                           [Mij; Mij],Ndof,Ndof);           
  end    
  %  Stiff = Stiff + sparse(double(elem2dof(:,i)),double(elem2dof(:,j)),Aij,Ndof,Ndof);    
        
    end
end
clear AA MM
%  
for p = 1:Quad2
     for i = 1:DG_base_num    
        for j = 1:DG_base_num
        Mij = (DG_phi(p,:)*DGcof(:,i))*(DG_phi(p,:)*DGcof(:,j))*ones(NT,1).*area;
        M_DG = M_DG + sparse(DG_dof(:,i),DG_dof(:,j),Mij,NDG,NDG);
        end
     end
end    
%% ç»„è£…åˆšåº¦çŸ©é˜µ 
%DPHI = zeros(10,10);

%Dxy = zeros(NT,base_num);
%tDxx = zeros(NT,3);
%tDyy = tDxx;
%tDxy = tDxx;

%testD = zeros(NT,10);
for p = 1:Quad2

Dx = zeros(NT,base_num);
Dy = zeros(NT,base_num);

l1 = lambda2(p,1);
l2 = lambda2(p,2);
l3 = lambda2(p,3);

for i = 1:base_num
    a = m7(i,1);
    b = m7(i,2);
    c = m7(i,3);

    if a > 0
        Dx(:,i) = Dx(:,i) + a * l1.^(a-1) .* l2.^b .* l3.^c * Dlambda(:,1,1);
        Dy(:,i) = Dy(:,i) + a * l1.^(a-1) .* l2.^b .* l3.^c * Dlambda(:,2,1);
    end
    if b > 0
        Dx(:,i) = Dx(:,i) + b * l1.^a .* l2.^(b-1) .* l3.^c * Dlambda(:,1,2);
        Dy(:,i) = Dy(:,i) + b * l1.^a .* l2.^(b-1) .* l3.^c * Dlambda(:,2,2);
    end
    if c > 0
        Dx(:,i) = Dx(:,i) + c * l1.^a .* l2.^b .* l3.^(c-1) * Dlambda(:,1,3);
        Dy(:,i) = Dy(:,i) + c * l1.^a .* l2.^b .* l3.^(c-1) * Dlambda(:,2,3);
    end
end


for i = 1:base_num
    for j = 1:DG_base_num
   Bij11  = weight2(p)*( dot(COF(i:base_num:end,:),Dx,2).*(DG_phi(p,:)*DGcof(:,j)*ones(NT,1))).*area;
   Bij22  = weight2(p)*( dot(COF(i:base_num:end,:),Dy,2).*(DG_phi(p,:)*DGcof(:,j)*ones(NT,1))).*area;
   Bij31  = weight2(p)*( dot(COF(i:base_num:end,:),Dy,2).*(DG_phi(p,:)*DGcof(:,j)*ones(NT,1))).*area;
   Bij32  = weight2(p)*( dot(COF(i:base_num:end,:),Dx,2).*(DG_phi(p,:)*DGcof(:,j)*ones(NT,1))).*area;
%  Bij11  = weight2(p)*( dot(COF(i:base_num:end,:),Dx,2).*(ones(NT,1))).*area;
%  Bij22  = weight2(p)*( dot(COF(i:base_num:end,:),Dy,2).*(ones(NT,1))).*area;


 BStiff = BStiff + sparse(double(elem2dof(:,i)),double(DG_dof(:,j)),Bij11,3*Ndof,2*NDG);
 BStiff = BStiff + sparse(Ndof + double(elem2dof(:,i)),NDG + double(DG_dof(:,j)),Bij22,3*Ndof,2*NDG);
 BStiff = BStiff + sparse(2*Ndof + double(elem2dof(:,i)),double(DG_dof(:,j)),Bij31,3*Ndof,2*NDG);
 BStiff = BStiff + sparse(2*Ndof + double(elem2dof(:,i)),NDG + double(DG_dof(:,j)),Bij32,3*Ndof,2*NDG); 

% BStiff = BStiff + sparse(double(elem2dof(:,i)),(1:NT)',Bij11,2*Ndof,NT);
% BStiff = BStiff + sparse(Ndof + double(elem2dof(:,i)),(1:NT)',Bij22,2*Ndof,NT);

 %BStiff = BStiff + sparse(double(elem2dof(:,i)),double(DG_dof(:,j)),Bij11,2*Ndof,NDG);
 %BStiff = BStiff + sparse(Ndof + double(elem2dof(:,i)),double(DG_dof(:,j)),Bij22,2*Ndof,NDG);
 %BStiff = BStiff + sparse(2*Ndof + double(elem2dof(:,i)),double(DG_dof(:,j)),Bij31,4*Ndof,2*NDG);
 %BStiff = BStiff + sparse(3*Ndof + double(elem2dof(:,i)),NDG + double(DG_dof(:,j)),Bij32,4*Ndof,2*NDG); 
     end
end


% 
% 
 end      
%Stiff = Stiff + (1/hsize^5)*pA;
%% å³ç«¯é¡?        
%big_Stiff = [iota*Stiff+M, sparse(Ndof,Ndof),sparse(Ndof,Ndof);sparse(Ndof,Ndof),iota*Stiff+M,sparse(Ndof,Ndof);sparse(Ndof,Ndof),sparse(Ndof,Ndof),2*iota*Stiff+2*M];%blkdiag(Stiff,Stiff);
%% ============================================================
% Build one proof-construction candidate sigma_proof
% sigma_proof = scalar_monomial * E on selected element(s)
%% ============================================================

coefProof = zeros(3*Ndof,1);

% Choose the element and monomial used in your proof construction.
% Example: lambda1^3 lambda2^2
t0 = 1;
expTarget = [3,2,0];

% Symmetric tensor direction E.
% Replace this by the tensor direction in your proof, e.g. t*t'
E = [1,0;0,0];

% find monomial index
idxMono = find(all(m5 == expTarget,2));
if isempty(idxMono)
    error('expTarget is not in m5. Check the polynomial degree or exponent.');
end

% local coefficient vector in monomial basis
eMono = zeros(base_num,1);
eMono(idxMono) = 1;

% local Hermite basis coefficients
COFt = COF((t0-1)*base_num+1:t0*base_num,:);
cLocal = (COFt') \ eMono;

% scatter to global stress vector
for i = 1:base_num
    gdof = elem2dof(t0,i);

    coefProof(gdof)          = coefProof(gdof)          + E(1,1)*cLocal(i);
    coefProof(Ndof + gdof)   = coefProof(Ndof + gdof)   + E(2,2)*cLocal(i);
    coefProof(2*Ndof + gdof) = coefProof(2*Ndof + gdof) + E(1,2)*cLocal(i);
end



beta_m = (2*mu+La)/(4*mu*(La+mu));
beta_G = -La/(4*mu*(La+mu));
big_Stiff = [beta_m*(iota*Stiff+M),beta_G*(iota*Stiff+M) ,sparse(Ndof,Ndof);beta_G*(iota*Stiff+M),beta_m*(iota*Stiff+M),sparse(Ndof,Ndof);sparse(Ndof,Ndof),sparse(Ndof,Ndof),(iota*Stiff+M)/mu];
All_Stiff = [big_Stiff,BStiff;BStiff',sparse(2*NDG,2*NDG)];


S = BStiff' * (big_Stiff \ BStiff); 
MStiff = [M_DG,zeros(NDG,NDG);zeros(NDG,NDG),M_DG];
S = (S+S')/2; %  symmetric the above result
M_DG = (M_DG+M_DG')/2;
  lam = eig(full(S), full(MStiff));
  beta = sqrt(min(lam));


end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


