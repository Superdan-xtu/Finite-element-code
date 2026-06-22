 
%% Introduction
% Note! this code now is test the case iota=0, i.e., is a HR-linear
% elasticity problem. 
% Just click run, error_s represent the L2 error of sigma
%  error1 is the L2 error of displacement
%%  Mesh Setting
for mesh_s = 1:7
 hsize = 1/(2*mesh_s);    
[node1,elem1] = squaremesh([0,1,0,1],hsize); 
refinesolut1 = load('solutiota10e4new.mat');
Rsolut = refinesolut1.solut;
refinemesh1 = load('elemla100000iota1.mat');
Relem = refinemesh1.elem;
refinenode1 = load('nodela100000iota1.mat');
Rnode = refinenode1.node;
refinestresscof = load('COFla100000iota1.mat');
Rcof = refinestresscof.COF;
refineDGcof = load('DGCOFla100000iota1.mat');
RDGcof = refineDGcof.DG_COF;
La = 1.0e+05;
mu = 0.3;
iota = 1.0e-02;
[RDlambda,Rarea,RelemSign] = gradbasis(Rnode,Relem);
[Relem2edge,Redge,Relem2edgeSign,RedgeSign,Redge2elem,Rbdedge] = dofedgemy(Relem);
%[ elemt , elemn ] = elemtau_n( node , Relem , Relem2edgeSign ,Relem2edge ,Redge );
%[Relem2dof1,Relem2edge,Redge,RbdDof,RfreeDof] = dofP4(Relem); 
[Relem2dof1,Relem2edge,Redge,RbdDof,RfreeDof] = C2dof(Relem,Rnode); 
Relem2dof = double(Relem2dof1);



% RefineN=1;
%  for k = 1:RefineN
% [elem2edge1,edge1,~,~,~,~] = dofedgemy(elem1);
% [ node1, elem1] = refinemeshme( node1,elem1,elem2edge1,edge1 );
% 
%  end

base_num = 15;
inter_node2 = 0.4*node1(elem1(:,1),:)+0.2*node1(elem1(:,2),:)+0.4*node1(elem1(:,3),:);
inter_node3 = 0.4*node1(elem1(:,1),:)+0.4*node1(elem1(:,2),:)+0.2*node1(elem1(:,3),:);
inter_node1 = 0.2*node1(elem1(:,1),:)+0.4*node1(elem1(:,2),:)+0.4*node1(elem1(:,3),:);
N = size(node1,1);
NT = size(elem1,1);
node = [node1;inter_node1;inter_node2;inter_node3];
elem_inter1 = [(N+1:N+NT)',(N+NT+1:N+2*NT)',(N+2*NT+1:N+3*NT)'];
elem_e1 = [elem1(:,2),elem1(:,3),(N+1:N+NT)' ];
elem_e2 = [(N+NT+1:N+2*NT)',elem1(:,3),elem1(:,1)];
elem_e3 = [(N+2*NT+1:N+3*NT)',elem1(:,1),elem1(:,2)];
elem_v1 = [ (N + NT + 1: N+ 2*NT)',elem1(:,1), (N+2*NT+1:N+3*NT)' ];
elem_v2 = [ (N + 1:N+NT)', (N+2*NT+1:N+3*NT)', elem1(:,2) ];
elem_v3 = [ (N + 1:N+NT)', elem1(:,3), (N+NT+1:N+2*NT)'  ];

elem = [elem_inter1;elem_v1;elem_v2;elem_v3;elem_e1;elem_e2;elem_e3];
Mesh_G = size(Relem,1)/size(elem,1);

[Dlambda,area,elemSign] = gradbasis(node,elem);
[elem2edge,edge,elem2edgeSign,edgeSign,edge2elem,bdedge] = dofedgemy(elem);
[ elemt , elemn ] = elemtau_n( node , elem , elem2edgeSign ,elem2edge ,edge );
elem = double(elem);
[elem2dof1,elem2edge,edge,bdDof,freeDof] = dofP4(elem); 
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
order2 = 9;
[lambda2,weight2] = quadpts(order2);
Quad2= size(lambda2,1);
phi = zeros(Quad2,base_num);

DG_base_num = 10;
DG_phi = zeros(Quad2,DG_base_num); 
DG_phi(:,1)=lambda2(:,1).^3;
DG_phi(:,2)=lambda2(:,2).^3;
DG_phi(:,3)=lambda2(:,3).^3;
DG_phi(:,4)=lambda2(:,1).^2.*lambda2(:,2);
DG_phi(:,5)=lambda2(:,1).^2.*lambda2(:,3);
DG_phi(:,6)=lambda2(:,2).^2.*lambda2(:,3);
DG_phi(:,7)=lambda2(:,1).*lambda2(:,2).^2;
DG_phi(:,8)=lambda2(:,1).*lambda2(:,3).^2;
DG_phi(:,9)=lambda2(:,2).*lambda2(:,3).^2;
DG_phi(:,10) = lambda2(:,1).*lambda2(:,2).*lambda2(:,3);

N = size(node,1);
NT = size(elem,1);

DG_dof = [(1:NT)', (NT+1:2*NT)', (2*NT+1:3*NT)',(3*NT+1:4*NT)',(4*NT+1:5*NT)',(5*NT+1:6*NT)',(6*NT+1:7*NT)',(7*NT+1:8*NT)',(8*NT+1:9*NT)',(9*NT+1:10*NT)'];
NDG = max(max(DG_dof));
BStiff = sparse(3*Ndof,2*NDG);
DGbase = zeros(DG_base_num,DG_base_num);
% DGbase(1,1) = 1;
% DGbase(2,2) = 1;
% DGbase(3,3) = 1;
DGm = [3,0,0;0,3,0;0,0,3;2,1,0;2,0,1;0,2,1;1,2,0;1,0,2;0,1,2;1,1,1];
for i = 1:DG_base_num
    for j = 1:DG_base_num
    DGbase(j,i) = com_quad((DGm(j,:)+DGm(i,:))');
    end
end
DGcof = inv(DGbase);
%DGcof = DGcof'; %% ???????????????
phi(:,1)=lambda2(:,1).^3;
phi(:,2)=lambda2(:,2).^3;
phi(:,3)=lambda2(:,3).^3;
phi(:,4)=lambda2(:,1).^2.*lambda2(:,2);
phi(:,5)=lambda2(:,1).^2.*lambda2(:,3);
phi(:,6)=lambda2(:,2).^2.*lambda2(:,3);
phi(:,7)=lambda2(:,1).*lambda2(:,2).^2;
phi(:,8)=lambda2(:,1).*lambda2(:,3).^2;
phi(:,9)=lambda2(:,2).*lambda2(:,3).^2;
phi(:,10) = lambda2(:,1).^2.*lambda2(:,2).*lambda2(:,3);
phi(:,11) = lambda2(:,1).*lambda2(:,2).^2.*lambda2(:,3);
phi(:,12) = lambda2(:,1).*lambda2(:,2).*lambda2(:,3).^2;
phi(:,13) = lambda2(:,1).^2.*lambda2(:,2).^2;
phi(:,15) = lambda2(:,2).^2.*lambda2(:,3).^2;
phi(:,14) = lambda2(:,3).^2.*lambda2(:,1).^2;

m3 = [3,0,0;0,3,0;0,0,3;2,1,0;2,0,1;0,2,1;1,2,0;1,0,2;0,1,2;2,1,1;1,2,1;1,1,2;2,2,0;2,0,2;0,2,2];
BA = zeros(3,15);
for i = 1:base_num
    BA(1,i) = com_quad((m3(i,:)+[1,0,0])');
    BA(2,i) = com_quad((m3(i,:)+[0,1,0])');
    BA(3,i) = com_quad((m3(i,:)+[0,0,1])');
end

BAE1 = zeros(3,15);
BAE2 = BAE1;
BAE3 = BAE1;
for i = 1:15
    BAE1(1,i) = 0^m3(i,1)*com_quad(m3(i,[2,3])'+[2;0]);
    BAE1(2,i) = 0^m3(i,1)*com_quad(m3(i,[2,3])'+[1;1]);
    BAE1(3,i) = 0^m3(i,1)*com_quad(m3(i,[2,3])'+[0;2]);
    
    BAE2(1,i) = 0^m3(i,2)*com_quad(m3(i,[1,3])'+[0;2]);
    BAE2(2,i) = 0^m3(i,2)*com_quad(m3(i,[1,3])'+[1;1]);
    BAE2(3,i) = 0^m3(i,2)*com_quad(m3(i,[1,3])'+[2;0]);
    
    BAE3(1,i) = 0^m3(i,3)*com_quad(m3(i,[1,2])'+[2;0]);
    BAE3(2,i) = 0^m3(i,3)*com_quad(m3(i,[1,2])'+[1;1]);
    BAE3(3,i) = 0^m3(i,3)*com_quad(m3(i,[1,2])'+[0;2]);
end


COF = zeros(base_num*NT,base_num);



for i = 1:NT
A = zeros(base_num,base_num);
A(1,1)=1;
A(2,2)=1;
A(3,3)=1;

A(4:6,:) = BAE1;
A(7:9,:) = BAE2;
A(10:12,:) = BAE3;
A(13:15,:)=BA;

cof = inv(A);
COF((i-1)*base_num+1:i*base_num,:)=cof'; 

end

Dx = zeros(NT,base_num);
Dy = zeros(NT,base_num);
AA = zeros(NT,(base_num*(base_num+1)/2));
MM = AA;


elem2dof = double(elem2dof);
for p = 1:Quad2
Dx(:,1) = 3*lambda2(p,1)^2*Dlambda(:,1,1);
Dy(:,1) = 3*lambda2(p,1)^2*Dlambda(:,2,1);

Dx(:,2) = 3*lambda2(p,2)^2*Dlambda(:,1,2);
Dy(:,2) = 3*lambda2(p,2)^2*Dlambda(:,2,2);

Dx(:,3) = 3*lambda2(p,3)^2*Dlambda(:,1,3);
Dy(:,3) = 3*lambda2(p,3)^2*Dlambda(:,2,3);

Dx(:,4) = 2*lambda2(p,1)*lambda2(p,2)*Dlambda(:,1,1)+lambda2(p,1)^2*Dlambda(:,1,2);
Dy(:,4) = 2*lambda2(p,1)*lambda2(p,2)*Dlambda(:,2,1)+lambda2(p,1)^2*Dlambda(:,2,2);

Dx(:,5) = 2*lambda2(p,1)*lambda2(p,3)*Dlambda(:,1,1)+lambda2(p,1)^2*Dlambda(:,1,3);
Dy(:,5) = 2*lambda2(p,1)*lambda2(p,3)*Dlambda(:,2,1)+lambda2(p,1)^2*Dlambda(:,2,3);

Dx(:,6) = 2*lambda2(p,3)*lambda2(p,2)*Dlambda(:,1,2)+lambda2(p,2)^2*Dlambda(:,1,3);
Dy(:,6) = 2*lambda2(p,3)*lambda2(p,2)*Dlambda(:,2,2)+lambda2(p,2)^2*Dlambda(:,2,3);

Dx(:,7) = 2*lambda2(p,1)*lambda2(p,2)*Dlambda(:,1,2)+lambda2(p,2)^2*Dlambda(:,1,1);
Dy(:,7) = 2*lambda2(p,1)*lambda2(p,2)*Dlambda(:,2,2)+lambda2(p,2)^2*Dlambda(:,2,1);

Dx(:,8) = 2*lambda2(p,1)*lambda2(p,3)*Dlambda(:,1,3)+lambda2(p,3)^2*Dlambda(:,1,1);
Dy(:,8) = 2*lambda2(p,1)*lambda2(p,3)*Dlambda(:,2,3)+lambda2(p,3)^2*Dlambda(:,2,1);

Dx(:,9) = 2*lambda2(p,3)*lambda2(p,2)*Dlambda(:,1,3)+lambda2(p,3)^2*Dlambda(:,1,2);
Dy(:,9) = 2*lambda2(p,3)*lambda2(p,2)*Dlambda(:,2,3)+lambda2(p,3)^2*Dlambda(:,2,2);
L = [1,2,3];
Dx(:,10) = 2*lambda2(p,L(1))*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,1,L(1)) + lambda2(p,L(1))^2*lambda2(p,L(3))*Dlambda(:,1,L(2))....
         + lambda2(p,L(1))^2*lambda2(p,L(2))*Dlambda(:,1,L(3));
Dy(:,10) = 2*lambda2(p,L(1))*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,2,L(1)) + lambda2(p,L(1))^2*lambda2(p,L(3))*Dlambda(:,2,L(2))....
         + lambda2(p,L(1))^2*lambda2(p,L(2))*Dlambda(:,2,L(3));
L = [2,1,3];
Dx(:,11) = 2*lambda2(p,L(1))*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,1,L(1)) + lambda2(p,L(1))^2*lambda2(p,L(3))*Dlambda(:,1,L(2))....
         + lambda2(p,L(1))^2*lambda2(p,L(2))*Dlambda(:,1,L(3));
Dy(:,11) = 2*lambda2(p,L(1))*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,2,L(1)) + lambda2(p,L(1))^2*lambda2(p,L(3))*Dlambda(:,2,L(2))....
         + lambda2(p,L(1))^2*lambda2(p,L(2))*Dlambda(:,2,L(3));    
L = [3,1,2];
Dx(:,12) = 2*lambda2(p,L(1))*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,1,L(1)) + lambda2(p,L(1))^2*lambda2(p,L(3))*Dlambda(:,1,L(2))....
         + lambda2(p,L(1))^2*lambda2(p,L(2))*Dlambda(:,1,L(3));
Dy(:,12) = 2*lambda2(p,L(1))*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,2,L(1)) + lambda2(p,L(1))^2*lambda2(p,L(3))*Dlambda(:,2,L(2))....
         + lambda2(p,L(1))^2*lambda2(p,L(2))*Dlambda(:,2,L(3));
Dx(:,13) = 2*lambda2(p,1)*lambda2(p,2)^2*Dlambda(:,1,1) + 2*lambda2(p,1)^2*lambda2(p,2)*Dlambda(:,1,2);
Dy(:,13) = 2*lambda2(p,1)*lambda2(p,2)^2*Dlambda(:,2,1) + 2*lambda2(p,1)^2*lambda2(p,2)*Dlambda(:,2,2);

Dx(:,15) = 2*lambda2(p,2)*lambda2(p,3)^2*Dlambda(:,1,2) + 2*lambda2(p,2)^2*lambda2(p,3)*Dlambda(:,1,3);
Dy(:,15) = 2*lambda2(p,2)*lambda2(p,3)^2*Dlambda(:,2,2) + 2*lambda2(p,2)^2*lambda2(p,3)*Dlambda(:,2,3);

Dx(:,14) = 2*lambda2(p,3)*lambda2(p,1)^2*Dlambda(:,1,3) + 2*lambda2(p,3)^2*lambda2(p,1)*Dlambda(:,1,1);
Dy(:,14) = 2*lambda2(p,3)*lambda2(p,1)^2*Dlambda(:,2,3) + 2*lambda2(p,3)^2*lambda2(p,1)*Dlambda(:,2,1);

    
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
   
%% ç»„è£…åˆšåº¦çŸ©é˜µ 
%DPHI = zeros(10,10);
Dx = zeros(NT,base_num);
Dy = zeros(NT,base_num);
%Dxy = zeros(NT,base_num);
%tDxx = zeros(NT,3);
%tDyy = tDxx;
%tDxy = tDxx;

%testD = zeros(NT,10);
for p = 1:Quad2

Dx(:,1) = 3*lambda2(p,1)^2*Dlambda(:,1,1);
Dy(:,1) = 3*lambda2(p,1)^2*Dlambda(:,2,1);

Dx(:,2) = 3*lambda2(p,2)^2*Dlambda(:,1,2);
Dy(:,2) = 3*lambda2(p,2)^2*Dlambda(:,2,2);

Dx(:,3) = 3*lambda2(p,3)^2*Dlambda(:,1,3);
Dy(:,3) = 3*lambda2(p,3)^2*Dlambda(:,2,3);

Dx(:,4) = 2*lambda2(p,1)*lambda2(p,2)*Dlambda(:,1,1)+lambda2(p,1)^2*Dlambda(:,1,2);
Dy(:,4) = 2*lambda2(p,1)*lambda2(p,2)*Dlambda(:,2,1)+lambda2(p,1)^2*Dlambda(:,2,2);

Dx(:,5) = 2*lambda2(p,1)*lambda2(p,3)*Dlambda(:,1,1)+lambda2(p,1)^2*Dlambda(:,1,3);
Dy(:,5) = 2*lambda2(p,1)*lambda2(p,3)*Dlambda(:,2,1)+lambda2(p,1)^2*Dlambda(:,2,3);

Dx(:,6) = 2*lambda2(p,3)*lambda2(p,2)*Dlambda(:,1,2)+lambda2(p,2)^2*Dlambda(:,1,3);
Dy(:,6) = 2*lambda2(p,3)*lambda2(p,2)*Dlambda(:,2,2)+lambda2(p,2)^2*Dlambda(:,2,3);

Dx(:,7) = 2*lambda2(p,1)*lambda2(p,2)*Dlambda(:,1,2)+lambda2(p,2)^2*Dlambda(:,1,1);
Dy(:,7) = 2*lambda2(p,1)*lambda2(p,2)*Dlambda(:,2,2)+lambda2(p,2)^2*Dlambda(:,2,1);

Dx(:,8) = 2*lambda2(p,1)*lambda2(p,3)*Dlambda(:,1,3)+lambda2(p,3)^2*Dlambda(:,1,1);
Dy(:,8) = 2*lambda2(p,1)*lambda2(p,3)*Dlambda(:,2,3)+lambda2(p,3)^2*Dlambda(:,2,1);

Dx(:,9) = 2*lambda2(p,3)*lambda2(p,2)*Dlambda(:,1,3)+lambda2(p,3)^2*Dlambda(:,1,2);
Dy(:,9) = 2*lambda2(p,3)*lambda2(p,2)*Dlambda(:,2,3)+lambda2(p,3)^2*Dlambda(:,2,2);
L = [1,2,3];
Dx(:,10) = 2*lambda2(p,L(1))*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,1,L(1)) + lambda2(p,L(1))^2*lambda2(p,L(3))*Dlambda(:,1,L(2))....
         + lambda2(p,L(1))^2*lambda2(p,L(2))*Dlambda(:,1,L(3));
Dy(:,10) = 2*lambda2(p,L(1))*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,2,L(1)) + lambda2(p,L(1))^2*lambda2(p,L(3))*Dlambda(:,2,L(2))....
         + lambda2(p,L(1))^2*lambda2(p,L(2))*Dlambda(:,2,L(3));
L = [2,1,3];
Dx(:,11) = 2*lambda2(p,L(1))*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,1,L(1)) + lambda2(p,L(1))^2*lambda2(p,L(3))*Dlambda(:,1,L(2))....
         + lambda2(p,L(1))^2*lambda2(p,L(2))*Dlambda(:,1,L(3));
Dy(:,11) = 2*lambda2(p,L(1))*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,2,L(1)) + lambda2(p,L(1))^2*lambda2(p,L(3))*Dlambda(:,2,L(2))....
         + lambda2(p,L(1))^2*lambda2(p,L(2))*Dlambda(:,2,L(3));    
L = [3,1,2];
Dx(:,12) = 2*lambda2(p,L(1))*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,1,L(1)) + lambda2(p,L(1))^2*lambda2(p,L(3))*Dlambda(:,1,L(2))....
         + lambda2(p,L(1))^2*lambda2(p,L(2))*Dlambda(:,1,L(3));
Dy(:,12) = 2*lambda2(p,L(1))*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,2,L(1)) + lambda2(p,L(1))^2*lambda2(p,L(3))*Dlambda(:,2,L(2))....
         + lambda2(p,L(1))^2*lambda2(p,L(2))*Dlambda(:,2,L(3));
Dx(:,13) = 2*lambda2(p,1)*lambda2(p,2)^2*Dlambda(:,1,1) + 2*lambda2(p,1)^2*lambda2(p,2)*Dlambda(:,1,2);
Dy(:,13) = 2*lambda2(p,1)*lambda2(p,2)^2*Dlambda(:,2,1) + 2*lambda2(p,1)^2*lambda2(p,2)*Dlambda(:,2,2);

Dx(:,15) = 2*lambda2(p,2)*lambda2(p,3)^2*Dlambda(:,1,2) + 2*lambda2(p,2)^2*lambda2(p,3)*Dlambda(:,1,3);
Dy(:,15) = 2*lambda2(p,2)*lambda2(p,3)^2*Dlambda(:,2,2) + 2*lambda2(p,2)^2*lambda2(p,3)*Dlambda(:,2,3);

Dx(:,14) = 2*lambda2(p,3)*lambda2(p,1)^2*Dlambda(:,1,3) + 2*lambda2(p,3)^2*lambda2(p,1)*Dlambda(:,1,1);
Dy(:,14) = 2*lambda2(p,3)*lambda2(p,1)^2*Dlambda(:,2,3) + 2*lambda2(p,3)^2*lambda2(p,1)*Dlambda(:,2,1);


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
beta_m = (2*mu+La)/(4*mu*(La+mu));
beta_G = -La/(4*mu*(La+mu));
big_Stiff = [beta_m*(iota^2*Stiff+M),beta_G*(iota^2*Stiff+M) ,sparse(Ndof,Ndof);beta_G*(iota^2*Stiff+M),beta_m*(iota^2*Stiff+M),sparse(Ndof,Ndof);sparse(Ndof,Ndof),sparse(Ndof,Ndof),(iota^2*Stiff+M)/mu];
All_Stiff = [big_Stiff,BStiff;BStiff',sparse(2*NDG,2*NDG)];

%% Right hand side
Bt = zeros(NT,DG_base_num);
Bt1 = Bt;
 for p  =1:Quad2        
  pxy = lambda2(p,1)*node(elem(:,1),:)+lambda2(p,2)*node(elem(:,2),:)+lambda2(p,3)*node(elem(:,3),:);   
  %[f1,f2]=Ipf_Straingradient(pxy,mu,La,ell);
  [f,f1] = La_Ipf_LSG(pxy,La,mu);
%[f,f1] = layer_Ipf_LSG(pxy,La,mu);
 % f = ones(NT,1);
% f1 = ones(NT,1);
 %   for i = 1:base_num
%   Bt(:,i)=Bt(:,i)+weight2(p)*(COF(i:base_num:end,:)*(phi(p,:)')).*f;  
%   Bt1(:,i)=Bt1(:,i)+weight2(p)*(COF(i:base_num:end,:)*(phi(p,:)')).*f1; 
%   end
  
for i = 1:DG_base_num
  Bt(:,i)=Bt(:,i)+weight2(p)*(DG_phi(p,:)*DGcof(:,i)*ones(NT,1)).*f;  
  Bt1(:,i)=Bt1(:,i)+weight2(p)*(DG_phi(p,:)*DGcof(:,i)*ones(NT,1)).*f1; 
end

 end
    Bt = Bt.*repmat(area,1,DG_base_num);
    Bt1 = Bt1.*repmat(area,1,DG_base_num);
    b = accumarray(DG_dof(:),Bt(:),[NDG 1]);
    b1 = accumarray(DG_dof(:),Bt1(:),[NDG 1]);
    right = [zeros(3*Ndof,1);b;b1];
  %  right = [zeros(3*Ndof,1);ones(2*NDG,1)];
    solut = All_Stiff\right;

%% error estimation
solut1 = solut(3*Ndof+1:3*Ndof+NDG);
solut_s= solut(1:Ndof);
solut_s1 = solut(Ndof+1:2*Ndof);
solut_s2 = solut(2*Ndof+1:3*Ndof);
solut2 = solut(3*Ndof+NDG+1:3*Ndof+2*NDG);
int_u = solut1(DG_dof);
int_u1 = solut2(DG_dof);
int_s = solut_s(elem2dof);
int_s1 = solut_s1(elem2dof);
int_s2 = solut_s2(elem2dof);
%intu = uh(elem2dof);
% TRnode = Rnode;
% TRelem = Relem;
% TCnode = node;
% TCelem = elem;
% TDlambda = Dlambda;
% TRDlambda = RDlambda;
% clear Relem Rnode RDlambda node elem Dlambda 
% Relem = TCelem;
% Rnode = TCnode;
% elem = TRelem;
% node = TRnode;
% %Dlamda = TRDlambda;
% %RDlambda = TDlambda;
%[RDlambda,Rarea,RelemSign] = gradbasis(Rnode,Relem);

RNT = size(Relem,1);
RNdof = max(max(Relem2dof));
%RDx = zeros(RNT,base_num);
%RDy = zeros(RNT,base_num);
RDG_base_num = 15;
Rbase_num = 21;
error = zeros(NT,1);
enor = zeros(RNT,1);
error_s = zeros(NT,1);
error_u = error_s;
Gerror_s = error_s;
Rsolut_s =  Rsolut(1:RNdof);
Rsolut_s1 = Rsolut(RNdof+1:2*RNdof);
Rsolut_s2 = Rsolut(2*RNdof+1:3*RNdof);
RN = size(Rnode,1);
RDG_dof = [Relem, RN+Relem, 2*RN+Relem, (3*RN+1:3*RN+RNT)', (3*RN+RNT+1:3*RN+2*RNT)',(3*RN+2*RNT+1:3*RN+3*RNT)',(3*RN+3*RNT+1:3*RN+4*RNT)',(3*RN+4*RNT+1:3*RN+5*RNT)',(3*RN+5*RNT+1:3*RN+6*RNT)'];
RNDG = max(max(RDG_dof));
NDG = max(max(DG_dof));
Rint_s = Rsolut_s(Relem2dof);
Rint_s1 = Rsolut_s1(Relem2dof);
Rint_s2 = Rsolut_s2(Relem2dof);
Rsolut_u = Rsolut(3*RNdof+1:3*RNdof+RNDG);
Rsolut_u1 = Rsolut(3*RNdof+RNDG+1:3*RNdof+2*RNDG);

%RDGdof = [(1:RNT)', (RNT+1:2*RNT)', (2*RNT+1:3*RNT)',(3*RNT+1:4*RNT)',(4*RNT+1:5*RNT)',(5*RNT+1:6*RNT)',(6*RNT+1:7*RNT)',(7*RNT+1:8*RNT)',(8*RNT+1:9*RNT)',(9*RNT+1:10*RNT)'];



Rint_u = Rsolut_u(RDG_dof);
Rint_u1 = Rsolut_u1(RDG_dof);

%% Gradient of Refine mesh



for p = 1:Quad2
    Ripinter = zeros(NT,RDG_base_num);
    Ripinter1 = zeros(NT,Rbase_num);
    DxRipinter = zeros(NT,Rbase_num);
    %DxRipinter1 = zeros(RNT,base_num);
    DyRipinter = zeros(NT,Rbase_num);
    %DyRipinter1 = zeros(RNT,base_num);
    Cipinter   = zeros(NT,DG_base_num);
    DxCipinter = zeros(NT,base_num);
    DyCipinter = zeros(NT,base_num);
    Cipinter1 = zeros(NT,base_num);
    
    Ripinteru1 = zeros(NT,RDG_base_num);
    Ripinters1 = zeros(NT,Rbase_num);
    DxRipinters1 = zeros(NT,Rbase_num);
    %DxRipinter1 = zeros(RNT,base_num);
    DyRipinters1 = zeros(NT,Rbase_num);
    %DyRipinter1 = zeros(RNT,base_num);
    Cipinteru1   = zeros(NT,DG_base_num);
    DxCipinters1 = zeros(NT,base_num);
    DyCipinters1 = zeros(NT,base_num);
    Cipinters1 = zeros(NT,base_num);
    
    Ripinters2 = zeros(NT,Rbase_num);
    DxRipinters2 = zeros(NT,Rbase_num);
    %DxRipinter1 = zeros(RNT,base_num);
    DyRipinters2 = zeros(NT,Rbase_num);
    %DyRipinter1 = zeros(RNT,base_num);
    DxCipinters2 = zeros(NT,base_num);
    DyCipinters2 = zeros(NT,base_num);
    Cipinters2 = zeros(NT,base_num);

%% solution on Cmesh 
Dx(:,1) = 3*lambda2(p,1)^2*Dlambda(:,1,1);
Dy(:,1) = 3*lambda2(p,1)^2*Dlambda(:,2,1);

Dx(:,2) = 3*lambda2(p,2)^2*Dlambda(:,1,2);
Dy(:,2) = 3*lambda2(p,2)^2*Dlambda(:,2,2);

Dx(:,3) = 3*lambda2(p,3)^2*Dlambda(:,1,3);
Dy(:,3) = 3*lambda2(p,3)^2*Dlambda(:,2,3);

Dx(:,4) = 2*lambda2(p,1)*lambda2(p,2)*Dlambda(:,1,1)+lambda2(p,1)^2*Dlambda(:,1,2);
Dy(:,4) = 2*lambda2(p,1)*lambda2(p,2)*Dlambda(:,2,1)+lambda2(p,1)^2*Dlambda(:,2,2);

Dx(:,5) = 2*lambda2(p,1)*lambda2(p,3)*Dlambda(:,1,1)+lambda2(p,1)^2*Dlambda(:,1,3);
Dy(:,5) = 2*lambda2(p,1)*lambda2(p,3)*Dlambda(:,2,1)+lambda2(p,1)^2*Dlambda(:,2,3);

Dx(:,6) = 2*lambda2(p,3)*lambda2(p,2)*Dlambda(:,1,2)+lambda2(p,2)^2*Dlambda(:,1,3);
Dy(:,6) = 2*lambda2(p,3)*lambda2(p,2)*Dlambda(:,2,2)+lambda2(p,2)^2*Dlambda(:,2,3);

Dx(:,7) = 2*lambda2(p,1)*lambda2(p,2)*Dlambda(:,1,2)+lambda2(p,2)^2*Dlambda(:,1,1);
Dy(:,7) = 2*lambda2(p,1)*lambda2(p,2)*Dlambda(:,2,2)+lambda2(p,2)^2*Dlambda(:,2,1);

Dx(:,8) = 2*lambda2(p,1)*lambda2(p,3)*Dlambda(:,1,3)+lambda2(p,3)^2*Dlambda(:,1,1);
Dy(:,8) = 2*lambda2(p,1)*lambda2(p,3)*Dlambda(:,2,3)+lambda2(p,3)^2*Dlambda(:,2,1);

Dx(:,9) = 2*lambda2(p,3)*lambda2(p,2)*Dlambda(:,1,3)+lambda2(p,3)^2*Dlambda(:,1,2);
Dy(:,9) = 2*lambda2(p,3)*lambda2(p,2)*Dlambda(:,2,3)+lambda2(p,3)^2*Dlambda(:,2,2);
L = [1,2,3];
Dx(:,10) = 2*lambda2(p,L(1))*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,1,L(1)) + lambda2(p,L(1))^2*lambda2(p,L(3))*Dlambda(:,1,L(2))....
         + lambda2(p,L(1))^2*lambda2(p,L(2))*Dlambda(:,1,L(3));
Dy(:,10) = 2*lambda2(p,L(1))*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,2,L(1)) + lambda2(p,L(1))^2*lambda2(p,L(3))*Dlambda(:,2,L(2))....
         + lambda2(p,L(1))^2*lambda2(p,L(2))*Dlambda(:,2,L(3));
L = [2,1,3];
Dx(:,11) = 2*lambda2(p,L(1))*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,1,L(1)) + lambda2(p,L(1))^2*lambda2(p,L(3))*Dlambda(:,1,L(2))....
         + lambda2(p,L(1))^2*lambda2(p,L(2))*Dlambda(:,1,L(3));
Dy(:,11) = 2*lambda2(p,L(1))*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,2,L(1)) + lambda2(p,L(1))^2*lambda2(p,L(3))*Dlambda(:,2,L(2))....
         + lambda2(p,L(1))^2*lambda2(p,L(2))*Dlambda(:,2,L(3));    
L = [3,1,2];
Dx(:,12) = 2*lambda2(p,L(1))*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,1,L(1)) + lambda2(p,L(1))^2*lambda2(p,L(3))*Dlambda(:,1,L(2))....
         + lambda2(p,L(1))^2*lambda2(p,L(2))*Dlambda(:,1,L(3));
Dy(:,12) = 2*lambda2(p,L(1))*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,2,L(1)) + lambda2(p,L(1))^2*lambda2(p,L(3))*Dlambda(:,2,L(2))....
         + lambda2(p,L(1))^2*lambda2(p,L(2))*Dlambda(:,2,L(3));
Dx(:,13) = 2*lambda2(p,1)*lambda2(p,2)^2*Dlambda(:,1,1) + 2*lambda2(p,1)^2*lambda2(p,2)*Dlambda(:,1,2);
Dy(:,13) = 2*lambda2(p,1)*lambda2(p,2)^2*Dlambda(:,2,1) + 2*lambda2(p,1)^2*lambda2(p,2)*Dlambda(:,2,2);

Dx(:,15) = 2*lambda2(p,2)*lambda2(p,3)^2*Dlambda(:,1,2) + 2*lambda2(p,2)^2*lambda2(p,3)*Dlambda(:,1,3);
Dy(:,15) = 2*lambda2(p,2)*lambda2(p,3)^2*Dlambda(:,2,2) + 2*lambda2(p,2)^2*lambda2(p,3)*Dlambda(:,2,3);

Dx(:,14) = 2*lambda2(p,3)*lambda2(p,1)^2*Dlambda(:,1,3) + 2*lambda2(p,3)^2*lambda2(p,1)*Dlambda(:,1,1);
Dy(:,14) = 2*lambda2(p,3)*lambda2(p,1)^2*Dlambda(:,2,3) + 2*lambda2(p,3)^2*lambda2(p,1)*Dlambda(:,2,1);   
    for i = 1:DG_base_num
      Cipinter(:,i)= Cipinter(:,i)+int_u(:,i).*((DGcof(i,:)*(DG_phi(p,:)')*ones(NT,1)));
      Cipinteru1(:,i)= Cipinteru1(:,i)+int_u1(:,i).*((DGcof(i,:)*(DG_phi(p,:)')*ones(NT,1)));
   end  
   for i = 1:base_num
      %% sigma11
       Cipinter1(:,i) = Cipinter1(:,i)+int_s(:,i).*((phi(p,:)*cof(:,i))*ones(NT,1));
      DxCipinter(:,i)= DxCipinter(:,i)+int_s(:,i).*(Dx*cof(:,i));
      DyCipinter(:,i)= DyCipinter(:,i)+int_s(:,i).*(Dy*cof(:,i));
      %%  sigma22
      Cipinters1(:,i) = Cipinters1(:,i)+int_s1(:,i).*((phi(p,:)*cof(:,i))*ones(NT,1));
      DxCipinters1(:,i)= DxCipinters1(:,i)+int_s1(:,i).*(Dx*cof(:,i));
      DyCipinters1(:,i)= DyCipinters1(:,i)+int_s1(:,i).*(Dy*cof(:,i));
      %% sigma12
      Cipinters2(:,i) = Cipinters2(:,i)+int_s2(:,i).*((phi(p,:)*cof(:,i))*ones(NT,1));
      DxCipinters2(:,i)= DxCipinters2(:,i)+int_s2(:,i).*(Dx*cof(:,i));
      DyCipinters2(:,i)= DyCipinters2(:,i)+int_s2(:,i).*(Dy*cof(:,i));
      
   end  
   
 for k = 1:NT
   pxy = lambda2(p,1)*node(elem(k,1),:)+lambda2(p,2)*node(elem(k,2),:)+lambda2(p,3)*node(elem(k,3),:);  
   pxy = repmat(pxy, [RNT,1]);
  [ phi1,phi2,phi3 ] = Cursefunction( Rnode, Relem, Rarea, pxy );
  CCelem =  find(abs(phi1+phi2+phi3-1)<1.0e-10);
   %CCelem = CCelem1(1); 
[ phi1,phi2,phi3 ] = Cursefunction( Rnode, Relem(CCelem,:), Rarea(CCelem), pxy(1,:) );
Cphi = zeros(1,Rbase_num);
CDGphi = zeros(1,RDG_base_num);
CDxphi = zeros(1,Rbase_num);
CDyphi = zeros(1,Rbase_num);

CDxphi(:,1) = 3*phi1^2*RDlambda(CCelem,1,1);
CDyphi(:,1) = 3*phi1^2*RDlambda(CCelem,2,1);

CDxphi(:,2) = 3*phi2^2*RDlambda(CCelem,1,2);
CDyphi(:,2) = 3*phi2^2*RDlambda(CCelem,2,2);

CDxphi(:,3) = 3*phi3^2*RDlambda(CCelem,1,3);
CDyphi(:,3) = 3*phi3^2*RDlambda(CCelem,2,3);

CDxphi(:,4) = 2*phi1*phi2*RDlambda(CCelem,1,1)+phi1^2*RDlambda(CCelem,1,2);
CDyphi(:,4) = 2*phi1*phi2*RDlambda(CCelem,2,1)+phi1^2*RDlambda(CCelem,2,2);

CDxphi(:,5) = 2*phi1*phi3*RDlambda(CCelem,1,1)+phi1^2*RDlambda(CCelem,1,3);
CDyphi(:,5) = 2*phi1*phi3*RDlambda(CCelem,2,1)+phi1^2*RDlambda(CCelem,2,3);

CDxphi(:,6) = 2*phi3*phi2*RDlambda(CCelem,1,2)+phi2^2*RDlambda(CCelem,1,3);
CDyphi(:,6) = 2*phi3*phi2*RDlambda(CCelem,2,2)+phi2^2*RDlambda(CCelem,2,3);

CDxphi(:,7) = 2*phi1*phi2*RDlambda(CCelem,1,2)+phi2^2*RDlambda(CCelem,1,1);
CDyphi(:,7) = 2*phi1*phi2*RDlambda(CCelem,2,2)+phi2^2*RDlambda(CCelem,2,1);

CDxphi(:,8) = 2*phi1*phi3*RDlambda(CCelem,1,3)+phi3^2*RDlambda(CCelem,1,1);
CDyphi(:,8) = 2*phi1*phi3*RDlambda(CCelem,2,3)+phi3^2*RDlambda(CCelem,2,1);

CDxphi(:,9) = 2*phi3*phi2*RDlambda(CCelem,1,3)+phi3^2*RDlambda(CCelem,1,2);
CDyphi(:,9) = 2*phi3*phi2*RDlambda(CCelem,2,3)+phi3^2*RDlambda(CCelem,2,2);

CDxphi(:,10)= 2*phi1*phi2*phi3*RDlambda(CCelem,1,1)+phi1^2*phi3*RDlambda(CCelem,1,2)+....
          phi1^2*phi2*RDlambda(CCelem,1,3);
CDyphi(:,10)= 2*phi1*phi2*phi3*RDlambda(CCelem,2,1)+phi1^2*phi3*RDlambda(CCelem,2,2)+....
          phi1^2*phi2*RDlambda(CCelem,2,3);

CDxphi(:,11)= 2*phi1*phi2*phi3*RDlambda(CCelem,1,2)+phi2^2*phi3*RDlambda(CCelem,1,1)+....
          phi2^2*phi1*RDlambda(CCelem,1,3);
CDyphi(:,11)= 2*phi1*phi2*phi3*RDlambda(CCelem,2,2)+phi2^2*phi3*RDlambda(CCelem,2,1)+....
          phi2^2*phi1*RDlambda(CCelem,2,3);

CDxphi(:,12)= 2*phi1*phi2*phi3*RDlambda(CCelem,1,3)+phi3^2*phi2*RDlambda(CCelem,1,1)+....
          phi3^2*phi1*RDlambda(CCelem,1,2);
CDyphi(:,12)= 2*phi1*phi2*phi3*RDlambda(CCelem,2,3)+phi3^2*phi2*RDlambda(CCelem,2,1)+....
          phi3^2*phi1*RDlambda(CCelem,2,2);

CDxphi(:,13)= 2*phi1*phi2^2*RDlambda(CCelem,1,1)+2*phi1^2*phi2*RDlambda(CCelem,1,2);     
CDyphi(:,13)= 2*phi1*phi2^2*RDlambda(CCelem,2,1)+2*phi1^2*phi2*RDlambda(CCelem,2,2);         
      
CDxphi(:,14)= 3*phi1^2*phi2^2*RDlambda(CCelem,1,1)+2*phi1^3*phi2*RDlambda(CCelem,1,2);     
CDyphi(:,14)= 3*phi1^2*phi2^2*RDlambda(CCelem,2,1)+2*phi1^3*phi2*RDlambda(CCelem,2,2); 

CDxphi(:,15)= 2*phi2*phi3^2*RDlambda(CCelem,1,2)+2*phi2^2*phi3*RDlambda(CCelem,1,3);     
CDyphi(:,15)= 2*phi2*phi3^2*RDlambda(CCelem,2,2)+2*phi2^2*phi3*RDlambda(CCelem,2,3);         
      
CDxphi(:,16)= 3*phi2^2*phi3^2*RDlambda(CCelem,1,2)+2*phi2^3*phi3*RDlambda(CCelem,1,3);     
CDyphi(:,16)= 3*phi2^2*phi3^2*RDlambda(CCelem,2,2)+2*phi2^3*phi3*RDlambda(CCelem,2,3);

CDxphi(:,17)= 2*phi3*phi1^2*RDlambda(CCelem,1,3)+2*phi3^2*phi1*RDlambda(CCelem,1,1);     
CDyphi(:,17)= 2*phi3*phi1^2*RDlambda(CCelem,2,3)+2*phi3^2*phi1*RDlambda(CCelem,2,1);         
      
CDxphi(:,18)= 3*phi3^2*phi1^2*RDlambda(CCelem,1,3)+2*phi3^3*phi1*RDlambda(CCelem,1,1);     
CDyphi(:,18)= 3*phi3^2*phi1^2*RDlambda(CCelem,2,3)+2*phi3^3*phi1*RDlambda(CCelem,2,1);
L = [1,2,3];
CDxphi(:,19) = 2*phi1*phi2^2*phi3*RDlambda(CCelem,1,L(1))+2*phi1^2*phi2*phi3*RDlambda(CCelem,1,L(2))....
             + phi1^2*phi2^2*RDlambda(CCelem,1,L(3));
CDyphi(:,19) = 2*phi1*phi2^2*phi3*RDlambda(CCelem,2,L(1))+2*phi1^2*phi2*phi3*RDlambda(CCelem,2,L(2))....
             +  phi1^2*phi2^2*RDlambda(CCelem,2,L(3)); 
L = [3,1,2];
CDxphi(:,20) = 2*phi3*phi1^2*phi2*RDlambda(CCelem,1,L(1))+2*phi3^2*phi1*phi2*RDlambda(CCelem,1,L(2))....
             + phi3^2*phi1^2*RDlambda(CCelem,1,L(3));
CDyphi(:,20) = 2*phi3*phi1^2*phi2*RDlambda(CCelem,2,L(1))+2*phi3^2*phi1*phi2*RDlambda(CCelem,2,L(2))....
             +  phi3^2*phi1^2*RDlambda(CCelem,2,L(3));          
L = [2,3,1];
CDxphi(:,21) = 2*phi2*phi3^2*phi1*RDlambda(CCelem,1,L(1))+2*phi2^2*phi3*phi1*RDlambda(CCelem,1,L(2))....
             + phi2^2*phi3^2*RDlambda(CCelem,1,L(3));
CDyphi(:,21) = 2*phi2*phi3^2*phi1*RDlambda(CCelem,2,L(1))+2*phi2^2*phi3*phi1*RDlambda(CCelem,2,L(2))....
             +  phi2^2*phi3^2*RDlambda(CCelem,2,L(3));


Cphi(:,1) = phi1.^3;
Cphi(:,2) = phi2.^3;
Cphi(:,3) = phi3.^3;
Cphi(:,4) = phi1.^2.*phi2;
Cphi(:,5) = phi1.^2.*phi3;
Cphi(:,6) = phi2.^2.*phi3;
Cphi(:,7) = phi1.*phi2.^2;
Cphi(:,8) = phi1.*phi3.^2;
Cphi(:,9) = phi2.*phi3.^2;
Cphi(:,10) = phi1.^2.*phi2.*phi3;
Cphi(:,11) = phi1.*phi2.^2.*phi3;
Cphi(:,12) = phi1.*phi2.*phi3.^2;

Cphi(:,13) = phi1.^2.*phi2.^2;
Cphi(:,14) = phi1.^3.*phi2.^2;

Cphi(:,15) = phi2.^2.*phi3.^2;
Cphi(:,16) = phi2.^3.*phi3.^2;

Cphi(:,17) = phi3.^2.*phi1^2;
Cphi(:,18) = phi3.^3.*phi1^2;

Cphi(:,19) = phi1.^2.*phi2.^2*phi3;
Cphi(:,20) = phi1.^2.*phi2*phi3.^2;
Cphi(:,21) = phi1.*phi2.^2*phi3.^2;

CDGphi(:,1) = phi1.^3;
CDGphi(:,2) = phi2.^3;
CDGphi(:,3) = phi3.^3;
CDGphi(:,4)=phi1.^2.*phi2;
CDGphi(:,5)=phi1.^2.*phi3;
CDGphi(:,6)=phi2.^2.*phi3;
CDGphi(:,7)=phi1.*phi2.^2;
CDGphi(:,8)=phi1.*phi3.^2;
CDGphi(:,9)=phi2.*phi3.^2;
CDGphi(:,10) = phi1.^2.*phi2.*phi3;
CDGphi(:,11) = phi1.*phi2.^2.*phi3;
CDGphi(:,12) = phi1.*phi2.*phi3.^2;
CDGphi(:,13) = phi1.^2.*phi2.^2;
CDGphi(:,15) = phi2.^2.*phi3.^2;
CDGphi(:,14) = phi3.^2.*phi1.^2;



Temp_Rcof = Rcof((CCelem-1)*Rbase_num+1:CCelem*Rbase_num,:);
Temp_RDGcof = RDGcof((CCelem-1)*RDG_base_num+1:CCelem*RDG_base_num,:);
Temp_Rcof = Temp_Rcof';
Temp_RDGcof = Temp_RDGcof';
for i = 1:RDG_base_num    
Ripinter(k,i) = Ripinter(k,i) + Rint_u(CCelem,i).*(CDGphi*Temp_RDGcof(:,i));
Ripinteru1(k,i) = Ripinteru1(k,i) + Rint_u1(CCelem,i).*(CDGphi*Temp_RDGcof(:,i));
end

for i = 1:Rbase_num    
Ripinter1(k,i)  = Ripinter1(k,i)  + Rint_s(CCelem,i).*(Cphi*Temp_Rcof(:,i));
DxRipinter(k,i) = DxRipinter(k,i) + Rint_s(CCelem,i).*(CDxphi*Temp_Rcof(:,i));
DyRipinter(k,i) = DyRipinter(k,i) + Rint_s(CCelem,i).*(CDyphi*Temp_Rcof(:,i));
%% sigma 22
Ripinters1(k,i)  = Ripinters1(k,i)  + Rint_s1(CCelem,i).*(Cphi*Temp_Rcof(:,i));
DxRipinters1(k,i) = DxRipinters1(k,i) + Rint_s1(CCelem,i).*(CDxphi*Temp_Rcof(:,i));
DyRipinters1(k,i) = DyRipinters1(k,i) + Rint_s1(CCelem,i).*(CDyphi*Temp_Rcof(:,i));
%% sigma12
Ripinters2(k,i)  = Ripinters2(k,i)  + Rint_s2(CCelem,i).*(Cphi*Temp_Rcof(:,i));
DxRipinters2(k,i) = DxRipinters2(k,i) + Rint_s2(CCelem,i).*(CDxphi*Temp_Rcof(:,i));
DyRipinters2(k,i) = DyRipinters2(k,i) + Rint_s2(CCelem,i).*(CDyphi*Temp_Rcof(:,i));

end

end
   


  Ripi_s = sum(Ripinter1,2);
  Cipi_s = sum(Cipinter1,2);
  DxCipi_s = sum(DxCipinter,2);
  DyCipi_s = sum(DyCipinter,2);
  DxRipi_s = sum(DxRipinter,2);
  DyRipi_s = sum(DyRipinter,2);
  Ripi_s1 = sum(Ripinters1,2);
  Cipi_s1 = sum(Cipinters1,2);
  DxCipi_s1 = sum(DxCipinters1,2);
  DyCipi_s1 = sum(DyCipinters1,2);
  DxRipi_s1 = sum(DxRipinters1,2);
  DyRipi_s1 = sum(DyRipinters1,2);
  
  Ripi_s2 = sum(Ripinters2,2);
  Cipi_s2 = sum(Cipinters2,2);
  DxCipi_s2 = sum(DxCipinters2,2);
  DyCipi_s2 = sum(DyCipinters2,2);
  DxRipi_s2 = sum(DxRipinters2,2);
  DyRipi_s2 = sum(DyRipinters2,2);
  
  Ripi_u = sum( Ripinter,2);
  Cipi_u = sum( Cipinter,2);
  Ripi_u1 = sum( Ripinteru1,2);
  Cipi_u1 = sum( Cipinteru1,2);
  
%  enor = enor + weight2(p)*ipi.^2;
%  error = error + weight2(p)*(ru1-ipi).^2;%+weight2(p)*(ru2-ipi1).^2;
  error_s = error_s + weight2(p)*((Ripi_s-Cipi_s).^2+(Ripi_s1-Cipi_s1).^2+2*(Ripi_s2-Cipi_s2).^2);
%  Gerror_s = Gerror_s + weight2(p)*(iota*((DxRipi_s-DxCipi_s).^2 + (DyRipi_s-DyCipi_s).^2 + (DxRipi_s1-DxCipi_s1).^2 + (DyRipi_s1-DyCipi_s1).^2 + 2*(DxRipi_s2-DxCipi_s2).^2 + 2*(DyRipi_s2-DyCipi_s2).^2 )......
%        +   (DxRipi_s-DxCipi_s).^2 + (DyRipi_s1-DyCipi_s1).^2 + (DxRipi_s2-DxCipi_s2).^2 + (DyRipi_s2-DyCipi_s2) );    
Gerror_s = Gerror_s + weight2(p)*(iota^2*((DxRipi_s-DxCipi_s).^2 + (DyRipi_s-DyCipi_s).^2 + (DxRipi_s1-DxCipi_s1).^2 + (DyRipi_s1-DyCipi_s1).^2 + 2*(DxRipi_s2-DxCipi_s2).^2 + 2*(DyRipi_s2-DyCipi_s2).^2 )......
        +   (DxRipi_s+DyRipi_s2-DxCipi_s-DyCipi_s2).^2 + (DyRipi_s1+DxRipi_s2-DyCipi_s1-DxCipi_s2).^2  );     
%Gerror_s = Gerror_s + weight2(p)*iota*(DxCipi_s).^2;
error_u = error_u + weight2(p)*((Ripi_u-Cipi_u).^2 + (Ripi_u1-Cipi_u1).^2);
%   error = error + weight2(p)*(ipi).^2;
%   error_s = error_s + weight2(p)*(ipi_s).^2;
end

%error = error.*Rarea;
%enor = enor.*Rarea;
%enor = sqrt(sum(enor));
%error1 = sum(error);
%error1 = sqrt(error1);
error_s=error_s.*area;
errors = sum(error_s);
errors = sqrt(errors);
Gerror_s=Gerror_s.*area;
Gerrors = sum(Gerror_s);
Gerrors = sqrt(Gerrors)+errors;
error_u = error_u.*area;
erroru = sum(error_u);
erroru = sqrt(erroru);


disp(['Hsize:=' num2str(hsize)]);
disp(['L2 Error of u:=' num2str(erroru,'%e')]);
%disp(['L2 Error of sigma_11:=' num2str(errors,'%e')]);
disp(['H1 Error of sigma:=' num2str(Gerrors,'%e')]);

clear
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


