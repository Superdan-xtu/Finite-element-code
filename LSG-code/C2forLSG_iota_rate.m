clear   
%% Introduction
% This script tests the asymptotic behavior with respect to the small
% parameter iota for the linear stress gradient elasticity problem.
% The script can be run directly. The value of iota is set below.
%
% The boundary-layer case uses the layer_* functions.
% The compatible non-boundary-layer case uses the La_* functions.
% The computed energy norm error is printed in the MATLAB command window.
%%  Mesh Setting
for HH = 1:1
hsize = 1/64;
 [node,elem] = squaremesh([0,1,0,1],hsize); 
% RefineN=4;
%  for k = 1:RefineN
% [elem2edge1,edge1,~,~,~,~] = dofedgemy(elem1);
% [ node1, elem1] = refinemeshme( node1,elem1,elem2edge1,edge1 );
% 
%  end

iota = 1.0e-03;
La = 1;
mu = 0.3;
base_num = 21;
%inter_node2 = 0.4*node1(elem1(:,1),:)+0.2*node1(elem1(:,2),:)+0.4*node1(elem1(:,3),:);
%inter_node3 = 0.4*node1(elem1(:,1),:)+0.4*node1(elem1(:,2),:)+0.2*node1(elem1(:,3),:);
%inter_node1 = 0.2*node1(elem1(:,1),:)+0.4*node1(elem1(:,2),:)+0.4*node1(elem1(:,3),:);
N = size(node,1);
NT = size(elem,1);
% node = [node1;inter_node1;inter_node2;inter_node3];
% elem_inter1 = [(N+1:N+NT)',(N+NT+1:N+2*NT)',(N+2*NT+1:N+3*NT)'];
% elem_e1 = [elem1(:,2),elem1(:,3),(N+1:N+NT)' ];
% elem_e2 = [(N+NT+1:N+2*NT)',elem1(:,3),elem1(:,1)];
% elem_e3 = [(N+2*NT+1:N+3*NT)',elem1(:,1),elem1(:,2)];
% elem_v1 = [ (N + NT + 1: N+ 2*NT)',elem1(:,1), (N+2*NT+1:N+3*NT)' ];
% elem_v2 = [ (N + 1:N+NT)', (N+2*NT+1:N+3*NT)', elem1(:,2) ];
% elem_v3 = [ (N + 1:N+NT)', elem1(:,3), (N+NT+1:N+2*NT)'  ];

%elem = [elem_inter1;elem_v1;elem_v2;elem_v3;elem_e1;elem_e2;elem_e3];


[Dlambda,area,elemSign] = gradbasis(node,elem);
[elem2edge,edge,elem2edgeSign,edgeSign,edge2elem,bdedge] = dofedgemy(elem);
[ elemt , elemn ] = elemtau_n( node , elem , elem2edgeSign ,elem2edge ,edge );
elem = double(elem);
[elem2dof1,elem2edge,edge,bdDof,freeDof] = C2dof(elem,node); 
Ndof = double(max(elem2dof1(:)));
elem2dof = [elem2dof1];
Stiff = sparse(Ndof,Ndof);
M =  Stiff;
% 缁熶竴娉曞悜
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
tDx = zeros(NT,base_num,6);    
tDy = zeros(NT,base_num,6);        
tDxx = zeros(NT,base_num,6);        
tDyy = zeros(NT,base_num,6);
tDxy = zeros(NT,base_num,6);    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tlambda2 = [1,0,0,0,0.5,0.5;0,1,0,0.5,0,0.5;0,0,1,0.5,0.5,0]; 
tlambda2 = tlambda2';
for p = 1:6
tDxx(:,1,p)=6*tlambda2(p,1).*(Dlambda(:,1,1).^2);
%Dphi(2,1,:)
tDxx(:,2,p)=6*tlambda2(p,2).*(Dlambda(:,1,2).^2);
%Dphi(3,1,:)
tDxx(:,3,p)=6*tlambda2(p,3).*(Dlambda(:,1,3).^2);
%Dphi(1,2,:)
tDyy(:,1,p)=6*tlambda2(p,1).*(Dlambda(:,2,1).^2);
%Dphi(2,2,:)
tDyy(:,2,p)=6*tlambda2(p,2).*(Dlambda(:,2,2).^2);
%Dphi(3,2,:)
tDyy(:,3,p)=6*tlambda2(p,3).*(Dlambda(:,2,3).^2);
tDxy(:,1,p)=6*tlambda2(p,1).*(Dlambda(:,1,1).*Dlambda(:,2,1));
tDxy(:,2,p)=6*tlambda2(p,2).*(Dlambda(:,1,2).*Dlambda(:,2,2));
tDxy(:,3,p)=6*tlambda2(p,3).*(Dlambda(:,1,3).*Dlambda(:,2,3));


%Dphi(4,1,:)
tDxx(:,4,p)=   2*tlambda2(p,2).*(Dlambda(:,1,1).^2)+4*tlambda2(p,1).*Dlambda(:,1,1).*Dlambda(:,1,2);
%Dphi(4,2,:)
tDyy(:,4,p)=   2*tlambda2(p,2).*(Dlambda(:,2,1).^2)+4*tlambda2(p,1).*Dlambda(:,2,1).*Dlambda(:,2,2);
%Dphi(4,3,:)
tDxy(:,4,p)=  2*tlambda2(p,2).*(Dlambda(:,1,1).*Dlambda(:,2,1))+....
            2*tlambda2(p,1).*(Dlambda(:,1,1).*Dlambda(:,2,2)+Dlambda(:,1,2).*Dlambda(:,2,1));

%Dphi(5,1,:)
tDxx(:,5,p) = 2*tlambda2(p,3).*(Dlambda(:,1,1).^2)+4*tlambda2(p,1).*Dlambda(:,1,1).*Dlambda(:,1,3);
%Dphi(5,2,:)
tDyy(:,5,p) = 2*tlambda2(p,3).*(Dlambda(:,2,1).^2)+4*tlambda2(p,1).*Dlambda(:,2,1).*Dlambda(:,2,3);
%Dphi(5,3,:)
tDxy(:,5,p)=  2*tlambda2(p,3).*(Dlambda(:,1,1).*Dlambda(:,2,1))+....
            2*tlambda2(p,1).*(Dlambda(:,1,1).*Dlambda(:,2,3)+Dlambda(:,1,3).*Dlambda(:,2,1));

%Dphi(6,1,:)
tDxx(:,6,p)=   2*tlambda2(p,3).*(Dlambda(:,1,2).^2)+4*tlambda2(p,2).*Dlambda(:,1,2).*Dlambda(:,1,3);
%Dphi(6,2,:)
tDyy(:,6,p)=   2*tlambda2(p,3).*(Dlambda(:,2,2).^2)+4*tlambda2(p,2).*Dlambda(:,2,2).*Dlambda(:,2,3);
%Dphi(6,3,:)
tDxy(:,6,p)=  2*tlambda2(p,3).*(Dlambda(:,1,2).*Dlambda(:,2,2))+....
            2*tlambda2(p,2).*(Dlambda(:,1,2).*Dlambda(:,2,3)+Dlambda(:,1,3).*Dlambda(:,2,2));

%Dphi(7,1,:)
tDxx(:,7,p)= 2*tlambda2(p,1).*(Dlambda(:,1,2).^2)+4*tlambda2(p,2).*Dlambda(:,1,2).*Dlambda(:,1,1);
%Dphi(7,2,:)
tDyy(:,7,p)= 2*tlambda2(p,1).*(Dlambda(:,2,2).^2)+4*tlambda2(p,2).*Dlambda(:,2,2).*Dlambda(:,2,1);
%Dphi(7,3,:)
tDxy(:,7,p)=  2*tlambda2(p,1).*(Dlambda(:,1,2).*Dlambda(:,2,2))+....
            2*tlambda2(p,2).*(Dlambda(:,1,2).*Dlambda(:,2,1)+Dlambda(:,1,1).*Dlambda(:,2,2));

%Dphi(8,1,:)
tDxx(:,8,p)=   2*tlambda2(p,1).*(Dlambda(:,1,3).^2)+4*tlambda2(p,3).*Dlambda(:,1,3).*Dlambda(:,1,1);
%Dphi(8,2,:)
tDyy(:,8,p)=   2*tlambda2(p,1).*(Dlambda(:,2,3).^2)+4*tlambda2(p,3).*Dlambda(:,2,3).*Dlambda(:,2,1);
%Dphi(8,3,:)
tDxy(:,8,p)=  2*tlambda2(p,1).*(Dlambda(:,1,3).*Dlambda(:,2,3))+....
            2*tlambda2(p,3).*(Dlambda(:,1,3).*Dlambda(:,2,1)+Dlambda(:,1,1).*Dlambda(:,2,3));

%Dphi(9,1,:)
tDxx(:,9,p)=2*tlambda2(p,2).*(Dlambda(:,1,3).^2)+4*tlambda2(p,3).*Dlambda(:,1,3).*Dlambda(:,1,2);
%Dphi(9,2,:)
tDyy(:,9,p)=2*tlambda2(p,2).*(Dlambda(:,2,3).^2)+4*tlambda2(p,3).*Dlambda(:,2,3).*Dlambda(:,2,2);
%Dphi(9,3,:)
tDxy(:,9,p)=  2*tlambda2(p,2).*(Dlambda(:,1,3).*Dlambda(:,2,3))+....
            2*tlambda2(p,3).*(Dlambda(:,1,3).*Dlambda(:,2,2)+Dlambda(:,1,2).*Dlambda(:,2,3));

%Dphi(10,1,:)
tDxx(:,10,p) = 2*tlambda2(p,2)*tlambda2(p,3)*Dlambda(:,1,1).^2+2*tlambda2(p,1)^2*Dlambda(:,1,2).*Dlambda(:,1,3)+....
             4*tlambda2(p,1)*tlambda2(p,2)*Dlambda(:,1,1).*Dlambda(:,1,3)+ 4*tlambda2(p,1)*tlambda2(p,3)*Dlambda(:,1,1).*Dlambda(:,1,2);
%Dphi(10,2,:)
tDyy(:,10,p) = 2*tlambda2(p,2)*tlambda2(p,3)*Dlambda(:,2,1).^2+2*tlambda2(p,1)^2*Dlambda(:,2,2).*Dlambda(:,2,3)+....
             4*tlambda2(p,1)*tlambda2(p,2)*Dlambda(:,2,1).*Dlambda(:,2,3)+ 4*tlambda2(p,1)*tlambda2(p,3)*Dlambda(:,2,1).*Dlambda(:,2,2);      
%Dphi(10,3,:)
tDxy(:,10,p)=  2*tlambda2(p,2)*tlambda2(p,3)*Dlambda(:,1,1).*Dlambda(:,2,1)+tlambda2(p,1)^2*(Dlambda(:,1,2).*Dlambda(:,2,3)+Dlambda(:,2,2).*Dlambda(:,1,3))+....
            2*tlambda2(p,1)*tlambda2(p,3)*(Dlambda(:,1,1).*Dlambda(:,2,2)+Dlambda(:,2,1).*Dlambda(:,1,2))+....
            2*tlambda2(p,1)*tlambda2(p,2)*(Dlambda(:,1,1).*Dlambda(:,2,3)+Dlambda(:,2,1).*Dlambda(:,1,3));
        
        %Dphi(10,1,:)
tDxx(:,11,p) =2*tlambda2(p,1)*tlambda2(p,3)*Dlambda(:,1,2).^2+2*tlambda2(p,2)^2*Dlambda(:,1,1).*Dlambda(:,1,3)+....
            4*tlambda2(p,1)*tlambda2(p,2)*Dlambda(:,1,2).*Dlambda(:,1,3)+ 4*tlambda2(p,2)*tlambda2(p,3)*Dlambda(:,1,1).*Dlambda(:,1,2);
%Dphi(10,2,:)
tDyy(:,11,p) = 2*tlambda2(p,1)*tlambda2(p,3)*Dlambda(:,2,2).^2+2*tlambda2(p,2)^2*Dlambda(:,2,1).*Dlambda(:,2,3)+....
             4*tlambda2(p,1)*tlambda2(p,2)*Dlambda(:,2,2).*Dlambda(:,2,3)+ 4*tlambda2(p,2)*tlambda2(p,3)*Dlambda(:,2,1).*Dlambda(:,2,2);      
%Dphi(10,3,:)
tDxy(:,11,p)= 2*tlambda2(p,1)*tlambda2(p,3)*Dlambda(:,1,2).*Dlambda(:,2,2)+tlambda2(p,2)^2*(Dlambda(:,1,1).*Dlambda(:,2,3)+Dlambda(:,2,1).*Dlambda(:,1,3))+....
            2*tlambda2(p,2)*tlambda2(p,3)*(Dlambda(:,1,1).*Dlambda(:,2,2)+Dlambda(:,2,1).*Dlambda(:,1,2))+....
            2*tlambda2(p,1)*tlambda2(p,2)*(Dlambda(:,1,2).*Dlambda(:,2,3)+Dlambda(:,2,2).*Dlambda(:,1,3));
        
%Dphi(10,1,:)
tDxx(:,12,p) = 2*tlambda2(p,2)*tlambda2(p,1)*Dlambda(:,1,3).^2+2*tlambda2(p,3)^2*Dlambda(:,1,2).*Dlambda(:,1,1)+....
             4*tlambda2(p,3)*tlambda2(p,2)*Dlambda(:,1,3).*Dlambda(:,1,1)+ 4*tlambda2(p,1)*tlambda2(p,3)*Dlambda(:,1,3).*Dlambda(:,1,2);
%Dphi(10,2,:)
tDyy(:,12,p) = 2*tlambda2(p,2)*tlambda2(p,1)*Dlambda(:,2,3).^2+2*tlambda2(p,3)^2*Dlambda(:,2,2).*Dlambda(:,2,1)+....
             4*tlambda2(p,3)*tlambda2(p,2)*Dlambda(:,2,1).*Dlambda(:,2,3)+ 4*tlambda2(p,1)*tlambda2(p,3)*Dlambda(:,2,3).*Dlambda(:,2,2);      
%Dphi(10,3,:)
tDxy(:,12,p)= 2*tlambda2(p,2)*tlambda2(p,1)*Dlambda(:,1,3).*Dlambda(:,2,3)+tlambda2(p,3)^2*(Dlambda(:,1,2).*Dlambda(:,2,1)+Dlambda(:,2,2).*Dlambda(:,1,1))+....
            2*tlambda2(p,1)*tlambda2(p,3)*(Dlambda(:,1,3).*Dlambda(:,2,2)+Dlambda(:,2,3).*Dlambda(:,1,2))+....
            2*tlambda2(p,3)*tlambda2(p,2)*(Dlambda(:,1,1).*Dlambda(:,2,3)+Dlambda(:,2,1).*Dlambda(:,1,3));  
        
tDxx(:,13,p)= 2*tlambda2(p,2)^2*Dlambda(:,1,1).^2+8*tlambda2(p,1)*tlambda2(p,2)*Dlambda(:,1,1).*Dlambda(:,1,2)+....
            2*tlambda2(p,1)^2*Dlambda(:,1,2).^2;
tDyy(:,13,p)= 2*tlambda2(p,2)^2*Dlambda(:,2,1).^2+8*tlambda2(p,1)*tlambda2(p,2)*Dlambda(:,2,1).*Dlambda(:,2,2)+....
            2*tlambda2(p,1)^2*Dlambda(:,2,2).^2;       
tDxy(:,13,p)= 2*tlambda2(p,2)^2*Dlambda(:,1,1).*Dlambda(:,2,1)+2*tlambda2(p,1)^2*Dlambda(:,1,2).*Dlambda(:,2,2)+....
            4*tlambda2(p,1)*tlambda2(p,2)*(Dlambda(:,1,1).*Dlambda(:,2,2)+Dlambda(:,2,1).*Dlambda(:,1,2));
       
tDxx(:,14,p)= 6*tlambda2(p,1)*tlambda2(p,2)^2*Dlambda(:,1,1).^2+12*tlambda2(p,1)^2*tlambda2(p,2)*Dlambda(:,1,1).*Dlambda(:,1,2)+....
            2*tlambda2(p,1)^3*Dlambda(:,1,2).^2;
tDyy(:,14,p)= 6*tlambda2(p,1)*tlambda2(p,2)^2*Dlambda(:,2,1).^2+12*tlambda2(p,1)^2*tlambda2(p,2)*Dlambda(:,2,1).*Dlambda(:,2,2)+....
            2*tlambda2(p,1)^3*Dlambda(:,2,2).^2;        
tDxy(:,14,p)= 6*tlambda2(p,1)*tlambda2(p,2)^2*Dlambda(:,1,1).*Dlambda(:,2,1)+2*tlambda2(p,1)^3*Dlambda(:,1,2).*Dlambda(:,2,2)+....
            6*tlambda2(p,1)^2*tlambda2(p,2)*(Dlambda(:,1,1).*Dlambda(:,2,2)+Dlambda(:,2,1).*Dlambda(:,1,2));
                            
tDxx(:,15,p)= 2*tlambda2(p,3)^2*Dlambda(:,1,2).^2+8*tlambda2(p,2)*tlambda2(p,3)*Dlambda(:,1,2).*Dlambda(:,1,3)+....
            2*tlambda2(p,2)^2*Dlambda(:,1,3).^2;
tDyy(:,15,p)= 2*tlambda2(p,3)^2*Dlambda(:,2,2).^2+8*tlambda2(p,2)*tlambda2(p,3)*Dlambda(:,2,2).*Dlambda(:,2,3)+....
            2*tlambda2(p,2)^2*Dlambda(:,2,3).^2;       
tDxy(:,15,p)= 2*tlambda2(p,3)^2*Dlambda(:,1,2).*Dlambda(:,2,2)+2*tlambda2(p,2)^2*Dlambda(:,1,3).*Dlambda(:,2,3)+....
            4*tlambda2(p,2)*tlambda2(p,3)*(Dlambda(:,1,2).*Dlambda(:,2,3)+Dlambda(:,2,2).*Dlambda(:,1,3));       
        
tDxx(:,16,p)= 6*tlambda2(p,2)*tlambda2(p,3)^2*Dlambda(:,1,2).^2+12*tlambda2(p,2)^2*tlambda2(p,3)*Dlambda(:,1,2).*Dlambda(:,1,3)+....
            2*tlambda2(p,2)^3*Dlambda(:,1,3).^2;
tDyy(:,16,p)= 6*tlambda2(p,2)*tlambda2(p,3)^2*Dlambda(:,2,2).^2+12*tlambda2(p,2)^2*tlambda2(p,3)*Dlambda(:,2,2).*Dlambda(:,2,3)+....
            2*tlambda2(p,2)^3*Dlambda(:,2,3).^2;       
tDxy(:,16,p)= 6*tlambda2(p,2)*tlambda2(p,3)^2*Dlambda(:,1,2).*Dlambda(:,2,2)+2*tlambda2(p,2)^3*Dlambda(:,1,3).*Dlambda(:,2,3)+....
            6*tlambda2(p,2)^2*tlambda2(p,3)*(Dlambda(:,1,2).*Dlambda(:,2,3)+Dlambda(:,2,2).*Dlambda(:,1,3));         
       
tDxx(:,17,p)= 2*tlambda2(p,1)^2*Dlambda(:,1,3).^2+8*tlambda2(p,3)*tlambda2(p,1)*Dlambda(:,1,3).*Dlambda(:,1,1)+....
            2*tlambda2(p,3)^2*Dlambda(:,1,1).^2;
tDyy(:,17,p)= 2*tlambda2(p,1)^2*Dlambda(:,2,3).^2+8*tlambda2(p,3)*tlambda2(p,1)*Dlambda(:,2,3).*Dlambda(:,2,1)+....
            2*tlambda2(p,3)^2*Dlambda(:,2,1).^2;       
tDxy(:,17,p)= 2*tlambda2(p,1)^2*Dlambda(:,1,3).*Dlambda(:,2,3)+2*tlambda2(p,3)^2*Dlambda(:,1,1).*Dlambda(:,2,1)+....
            4*tlambda2(p,3)*tlambda2(p,1)*(Dlambda(:,1,3).*Dlambda(:,2,1)+Dlambda(:,2,3).*Dlambda(:,1,1));
       
tDxx(:,18,p)= 6*tlambda2(p,3)*tlambda2(p,1)^2*Dlambda(:,1,3).^2+12*tlambda2(p,3)^2*tlambda2(p,1)*Dlambda(:,1,3).*Dlambda(:,1,1)+....
            2*tlambda2(p,3)^3*Dlambda(:,1,1).^2;
tDyy(:,18,p)= 6*tlambda2(p,3)*tlambda2(p,1)^2*Dlambda(:,2,3).^2+12*tlambda2(p,3)^2*tlambda2(p,1)*Dlambda(:,2,3).*Dlambda(:,2,1)+....
            2*tlambda2(p,3)^3*Dlambda(:,2,1).^2;       
tDxy(:,18,p)= 6*tlambda2(p,3)*tlambda2(p,1)^2*Dlambda(:,1,3).*Dlambda(:,2,3)+2*tlambda2(p,3)^3*Dlambda(:,1,1).*Dlambda(:,2,1)+....
            6*tlambda2(p,3)^2*tlambda2(p,1)*(Dlambda(:,1,3).*Dlambda(:,2,1)+Dlambda(:,2,3).*Dlambda(:,1,1));   
       
tDx(:,1,p) = 3*tlambda2(p,1)^2*Dlambda(:,1,1);
tDy(:,1,p) = 3*tlambda2(p,1)^2*Dlambda(:,2,1);

tDx(:,2,p) = 3*tlambda2(p,2)^2*Dlambda(:,1,2);
tDy(:,2,p) = 3*tlambda2(p,2)^2*Dlambda(:,2,2);

tDx(:,3,p) = 3*tlambda2(p,3)^2*Dlambda(:,1,3);
tDy(:,3,p) = 3*tlambda2(p,3)^2*Dlambda(:,2,3);

tDx(:,4,p) = 2*tlambda2(p,1)*tlambda2(p,2)*Dlambda(:,1,1)+tlambda2(p,1)^2*Dlambda(:,1,2);
tDy(:,4,p) = 2*tlambda2(p,1)*tlambda2(p,2)*Dlambda(:,2,1)+tlambda2(p,1)^2*Dlambda(:,2,2);

tDx(:,5,p) = 2*tlambda2(p,1)*tlambda2(p,3)*Dlambda(:,1,1)+tlambda2(p,1)^2*Dlambda(:,1,3);
tDy(:,5,p) = 2*tlambda2(p,1)*tlambda2(p,3)*Dlambda(:,2,1)+tlambda2(p,1)^2*Dlambda(:,2,3);

tDx(:,6,p) = 2*tlambda2(p,3)*tlambda2(p,2)*Dlambda(:,1,2)+tlambda2(p,2)^2*Dlambda(:,1,3);
tDy(:,6,p) = 2*tlambda2(p,3)*tlambda2(p,2)*Dlambda(:,2,2)+tlambda2(p,2)^2*Dlambda(:,2,3);

tDx(:,7,p) = 2*tlambda2(p,1)*tlambda2(p,2)*Dlambda(:,1,2)+tlambda2(p,2)^2*Dlambda(:,1,1);
tDy(:,7,p) = 2*tlambda2(p,1)*tlambda2(p,2)*Dlambda(:,2,2)+tlambda2(p,2)^2*Dlambda(:,2,1);

tDx(:,8,p) = 2*tlambda2(p,1)*tlambda2(p,3)*Dlambda(:,1,3)+tlambda2(p,3)^2*Dlambda(:,1,1);
tDy(:,8,p) = 2*tlambda2(p,1)*tlambda2(p,3)*Dlambda(:,2,3)+tlambda2(p,3)^2*Dlambda(:,2,1);

tDx(:,9,p) = 2*tlambda2(p,3)*tlambda2(p,2)*Dlambda(:,1,3)+tlambda2(p,3)^2*Dlambda(:,1,2);
tDy(:,9,p) = 2*tlambda2(p,3)*tlambda2(p,2)*Dlambda(:,2,3)+tlambda2(p,3)^2*Dlambda(:,2,2);

tDx(:,10,p)= 2*tlambda2(p,1)*tlambda2(p,2)*tlambda2(p,3)*Dlambda(:,1,1)+tlambda2(p,1)^2*tlambda2(p,3)*Dlambda(:,1,2)+....
             tlambda2(p,1)^2*tlambda2(p,2)*Dlambda(:,1,3);
tDy(:,10,p)= 2*tlambda2(p,1)*tlambda2(p,2)*tlambda2(p,3)*Dlambda(:,2,1)+tlambda2(p,1)^2*tlambda2(p,3)*Dlambda(:,2,2)+....
             tlambda2(p,1)^2*tlambda2(p,2)*Dlambda(:,2,3);

tDx(:,11,p)= 2*tlambda2(p,2)*tlambda2(p,1)*tlambda2(p,3)*Dlambda(:,1,2)+tlambda2(p,2)^2*tlambda2(p,3)*Dlambda(:,1,1)+....
            tlambda2(p,2)^2*tlambda2(p,1)*Dlambda(:,1,3);
tDy(:,11,p)= 2*tlambda2(p,2)*tlambda2(p,1)*tlambda2(p,3)*Dlambda(:,2,2)+tlambda2(p,2)^2*tlambda2(p,3)*Dlambda(:,2,1)+....
          tlambda2(p,2)^2*tlambda2(p,1)*Dlambda(:,2,3);

tDx(:,12,p)= 2*tlambda2(p,3)*tlambda2(p,1)*tlambda2(p,2)*Dlambda(:,1,3)+tlambda2(p,3)^2*tlambda2(p,2)*Dlambda(:,1,1)+....
             tlambda2(p,3)^2*tlambda2(p,1)*Dlambda(:,1,2);
tDy(:,12,p)= 2*tlambda2(p,3)*tlambda2(p,1)*tlambda2(p,2)*Dlambda(:,2,3)+tlambda2(p,3)^2*tlambda2(p,2)*Dlambda(:,2,1)+....
             tlambda2(p,3)^2*tlambda2(p,1)*Dlambda(:,2,2);

tDx(:,13,p)= 2*tlambda2(p,1)*tlambda2(p,2)^2*Dlambda(:,1,1)+2*tlambda2(p,1)^2*tlambda2(p,2)*Dlambda(:,1,2);     
tDy(:,13,p)= 2*tlambda2(p,1)*tlambda2(p,2)^2*Dlambda(:,2,1)+2*tlambda2(p,1)^2*tlambda2(p,2)*Dlambda(:,2,2);         
      
tDx(:,14,p)= 3*tlambda2(p,1)^2*tlambda2(p,2)^2*Dlambda(:,1,1)+2*tlambda2(p,1)^3*tlambda2(p,2)*Dlambda(:,1,2);     
tDy(:,14,p)= 3*tlambda2(p,1)^2*tlambda2(p,2)^2*Dlambda(:,2,1)+2*tlambda2(p,1)^3*tlambda2(p,2)*Dlambda(:,2,2); 

tDx(:,15,p)= 2*tlambda2(p,2)*tlambda2(p,3)^2*Dlambda(:,1,2)+2*tlambda2(p,2)^2*tlambda2(p,3)*Dlambda(:,1,3);     
tDy(:,15,p)= 2*tlambda2(p,2)*tlambda2(p,3)^2*Dlambda(:,2,2)+2*tlambda2(p,2)^2*tlambda2(p,3)*Dlambda(:,2,3);         
      
tDx(:,16,p)= 3*tlambda2(p,2)^2*tlambda2(p,3)^2*Dlambda(:,1,2)+2*tlambda2(p,2)^3*tlambda2(p,3)*Dlambda(:,1,3);     
tDy(:,16,p)= 3*tlambda2(p,2)^2*tlambda2(p,3)^2*Dlambda(:,2,2)+2*tlambda2(p,2)^3*tlambda2(p,3)*Dlambda(:,2,3);

tDx(:,17,p)= 2*tlambda2(p,3)*tlambda2(p,1)^2*Dlambda(:,1,3)+2*tlambda2(p,3)^2*tlambda2(p,1)*Dlambda(:,1,1);     
tDy(:,17,p)= 2*tlambda2(p,3)*tlambda2(p,1)^2*Dlambda(:,2,3)+2*tlambda2(p,3)^2*tlambda2(p,1)*Dlambda(:,2,1);         
      
tDx(:,18,p)= 3*tlambda2(p,3)^2*tlambda2(p,1)^2*Dlambda(:,1,3)+2*tlambda2(p,3)^3*tlambda2(p,1)*Dlambda(:,1,1);     
tDy(:,18,p)= 3*tlambda2(p,3)^2*tlambda2(p,1)^2*Dlambda(:,2,3)+2*tlambda2(p,3)^3*tlambda2(p,1)*Dlambda(:,2,1);       

L = [1,2,3];
tDx(:,19,p) = 2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*tlambda2(p,L(3))*Dlambda(:,1,L(1))+2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,1,L(2))....
             +  tlambda2(p,L(1))^2*tlambda2(p,L(2))^2*Dlambda(:,1,L(3));
tDy(:,19,p) = 2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*tlambda2(p,L(3))*Dlambda(:,2,L(1))+2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,2,L(2))....
             +  tlambda2(p,L(1))^2*tlambda2(p,L(2))^2*Dlambda(:,2,L(3)); 
tDxx(:,19,p) = 2*tlambda2(p,L(2))^2*tlambda2(p,L(3))*Dlambda(:,1,L(1)).^2 + 4*tlambda2(p,L(1))*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,1,L(1)).*Dlambda(:,1,L(2))+2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*Dlambda(:,1,L(1)).*Dlambda(:,1,L(3))....
             + 4*tlambda2(p,L(1))*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,1,L(2)).*Dlambda(:,1,L(1)) + 2*tlambda2(p,L(1))^2*tlambda2(p,L(3))*Dlambda(:,1,L(2)).^2 + 2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*Dlambda(:,1,L(2)).*Dlambda(:,1,L(3))....
             + 2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*Dlambda(:,1,L(3)).*Dlambda(:,1,L(1)) + 2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*Dlambda(:,1,L(3)).*Dlambda(:,1,L(2));
tDyy(:,19,p) = 2*tlambda2(p,L(2))^2*tlambda2(p,L(3))*Dlambda(:,2,L(1)).^2 + 4*tlambda2(p,L(1))*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,2,L(1)).*Dlambda(:,2,L(2))+2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*Dlambda(:,2,L(1)).*Dlambda(:,2,L(3))....
             + 4*tlambda2(p,L(1))*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,2,L(2)).*Dlambda(:,2,L(1)) + 2*tlambda2(p,L(1))^2*tlambda2(p,L(3))*Dlambda(:,2,L(2)).^2 + 2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*Dlambda(:,2,L(2)).*Dlambda(:,2,L(3))....
             + 2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*Dlambda(:,2,L(3)).*Dlambda(:,2,L(1)) + 2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*Dlambda(:,2,L(3)).*Dlambda(:,2,L(2));
tDxy(:,19,p) = 2*tlambda2(p,L(2))^2*tlambda2(p,L(3))*Dlambda(:,1,L(1)).*Dlambda(:,2,L(1)) + 4*tlambda2(p,L(1))*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,1,L(1)).*Dlambda(:,2,L(2)) + 2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*Dlambda(:,1,L(1)).*Dlambda(:,2,L(3))....
             + 4*tlambda2(p,L(1))*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,1,L(2)).*Dlambda(:,2,L(1)) +   2*tlambda2(p,L(1))^2*tlambda2(p,L(3))*Dlambda(:,1,L(2)).*Dlambda(:,2,L(2)) +  2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*Dlambda(:,1,L(2)).*Dlambda(:,2,L(3))....
             + 2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*Dlambda(:,1,L(3)).*Dlambda(:,2,L(1)) +   2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*Dlambda(:,1,L(3)).*Dlambda(:,2,L(2));
L = [3,1,2];
tDx(:,20,p) = 2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*tlambda2(p,L(3))*Dlambda(:,1,L(1))+2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,1,L(2))....
             +  tlambda2(p,L(1))^2*tlambda2(p,L(2))^2*Dlambda(:,1,L(3));
tDy(:,20,p) = 2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*tlambda2(p,L(3))*Dlambda(:,2,L(1))+2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,2,L(2))....
             +  tlambda2(p,L(1))^2*tlambda2(p,L(2))^2*Dlambda(:,2,L(3)); 
tDxx(:,20,p) = 2*tlambda2(p,L(2))^2*tlambda2(p,L(3))*Dlambda(:,1,L(1)).^2 + 4*tlambda2(p,L(1))*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,1,L(1)).*Dlambda(:,1,L(2))+2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*Dlambda(:,1,L(1)).*Dlambda(:,1,L(3))....
             + 4*tlambda2(p,L(1))*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,1,L(2)).*Dlambda(:,1,L(1)) + 2*tlambda2(p,L(1))^2*tlambda2(p,L(3))*Dlambda(:,1,L(2)).^2 + 2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*Dlambda(:,1,L(2)).*Dlambda(:,1,L(3))....
             + 2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*Dlambda(:,1,L(3)).*Dlambda(:,1,L(1)) + 2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*Dlambda(:,1,L(3)).*Dlambda(:,1,L(2));
tDyy(:,20,p) = 2*tlambda2(p,L(2))^2*tlambda2(p,L(3))*Dlambda(:,2,L(1)).^2 + 4*tlambda2(p,L(1))*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,2,L(1)).*Dlambda(:,2,L(2))+2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*Dlambda(:,2,L(1)).*Dlambda(:,2,L(3))....
             + 4*tlambda2(p,L(1))*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,2,L(2)).*Dlambda(:,2,L(1)) + 2*tlambda2(p,L(1))^2*tlambda2(p,L(3))*Dlambda(:,2,L(2)).^2 + 2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*Dlambda(:,2,L(2)).*Dlambda(:,2,L(3))....
             + 2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*Dlambda(:,2,L(3)).*Dlambda(:,2,L(1)) + 2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*Dlambda(:,2,L(3)).*Dlambda(:,2,L(2));
tDxy(:,20,p) = 2*tlambda2(p,L(2))^2*tlambda2(p,L(3))*Dlambda(:,1,L(1)).*Dlambda(:,2,L(1)) + 4*tlambda2(p,L(1))*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,1,L(1)).*Dlambda(:,2,L(2)) + 2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*Dlambda(:,1,L(1)).*Dlambda(:,2,L(3))....
             + 4*tlambda2(p,L(1))*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,1,L(2)).*Dlambda(:,2,L(1)) +   2*tlambda2(p,L(1))^2*tlambda2(p,L(3))*Dlambda(:,1,L(2)).*Dlambda(:,2,L(2)) +  2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*Dlambda(:,1,L(2)).*Dlambda(:,2,L(3))....
             + 2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*Dlambda(:,1,L(3)).*Dlambda(:,2,L(1)) +   2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*Dlambda(:,1,L(3)).*Dlambda(:,2,L(2));
L = [2,3,1];
tDx(:,21,p) = 2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*tlambda2(p,L(3))*Dlambda(:,1,L(1))+2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,1,L(2))....
             +  tlambda2(p,L(1))^2*tlambda2(p,L(2))^2*Dlambda(:,1,L(3));
tDy(:,21,p) = 2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*tlambda2(p,L(3))*Dlambda(:,2,L(1))+2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,2,L(2))....
             +  tlambda2(p,L(1))^2*tlambda2(p,L(2))^2*Dlambda(:,2,L(3)); 
tDxx(:,21,p) = 2*tlambda2(p,L(2))^2*tlambda2(p,L(3))*Dlambda(:,1,L(1)).^2 + 4*tlambda2(p,L(1))*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,1,L(1)).*Dlambda(:,1,L(2))+2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*Dlambda(:,1,L(1)).*Dlambda(:,1,L(3))....
             + 4*tlambda2(p,L(1))*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,1,L(2)).*Dlambda(:,1,L(1)) + 2*tlambda2(p,L(1))^2*tlambda2(p,L(3))*Dlambda(:,1,L(2)).^2 + 2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*Dlambda(:,1,L(2)).*Dlambda(:,1,L(3))....
             + 2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*Dlambda(:,1,L(3)).*Dlambda(:,1,L(1)) + 2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*Dlambda(:,1,L(3)).*Dlambda(:,1,L(2));
tDyy(:,21,p) = 2*tlambda2(p,L(2))^2*tlambda2(p,L(3))*Dlambda(:,2,L(1)).^2 + 4*tlambda2(p,L(1))*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,2,L(1)).*Dlambda(:,2,L(2))+2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*Dlambda(:,2,L(1)).*Dlambda(:,2,L(3))....
             + 4*tlambda2(p,L(1))*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,2,L(2)).*Dlambda(:,2,L(1)) + 2*tlambda2(p,L(1))^2*tlambda2(p,L(3))*Dlambda(:,2,L(2)).^2 + 2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*Dlambda(:,2,L(2)).*Dlambda(:,2,L(3))....
             + 2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*Dlambda(:,2,L(3)).*Dlambda(:,2,L(1)) + 2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*Dlambda(:,2,L(3)).*Dlambda(:,2,L(2));
tDxy(:,21,p) = 2*tlambda2(p,L(2))^2*tlambda2(p,L(3))*Dlambda(:,1,L(1)).*Dlambda(:,2,L(1)) + 4*tlambda2(p,L(1))*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,1,L(1)).*Dlambda(:,2,L(2)) + 2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*Dlambda(:,1,L(1)).*Dlambda(:,2,L(3))....
             + 4*tlambda2(p,L(1))*tlambda2(p,L(2))*tlambda2(p,L(3))*Dlambda(:,1,L(2)).*Dlambda(:,2,L(1)) +   2*tlambda2(p,L(1))^2*tlambda2(p,L(3))*Dlambda(:,1,L(2)).*Dlambda(:,2,L(2)) +  2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*Dlambda(:,1,L(2)).*Dlambda(:,2,L(3))....
             + 2*tlambda2(p,L(1))*tlambda2(p,L(2))^2*Dlambda(:,1,L(3)).*Dlambda(:,2,L(1)) +   2*tlambda2(p,L(1))^2*tlambda2(p,L(2))*Dlambda(:,1,L(3)).*Dlambda(:,2,L(2));
end       

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
phi(:,14) = lambda2(:,1).^3.*lambda2(:,2).^2;

phi(:,15) = lambda2(:,2).^2.*lambda2(:,3).^2;
phi(:,16) = lambda2(:,2).^3.*lambda2(:,3).^2;

phi(:,17) = lambda2(:,3).^2.*lambda2(:,1).^2;
phi(:,18) = lambda2(:,3).^3.*lambda2(:,1).^2;

phi(:,19) = lambda2(:,1).^2.*lambda2(:,2).^2.*lambda2(:,3);
phi(:,20) = lambda2(:,1).^2.*lambda2(:,2).*lambda2(:,3).^2;
phi(:,21) = lambda2(:,1).*lambda2(:,2).^2.*lambda2(:,3).^2;

m3 = [3,0,0;0,3,0;0,0,3;2,1,0;2,0,1;0,2,1;1,2,0;1,0,2;0,1,2;2,1,1;1,2,1;1,1,2;2,2,0;3,2,0;0,2,2;0,3,2;2,0,2;2,0,3;2,2,1;2,1,2;1,2,2];
COF = zeros(base_num*NT,base_num);
A = zeros(base_num,base_num);
A(1,1) = 1;
A(2,2) = 1;
A(3,3) = 1;
for i = 1:base_num
A(19,i) = com_quad(([1,0,0]+m3(i,:))'); 
A(20,i) = com_quad(([0,1,0]+m3(i,:))'); 
A(21,i) = com_quad(([0,0,1]+m3(i,:))'); 
end
for i = 1:NT
% nor = zeros(3,2);
% nor(1,:)=elemn(i,:,1);
% nor(2,:)=elemn(i,:,2);
% nor(3,:)=elemn(i,:,3);
A(4:18,:) = [tDx(i,:,1);tDx(i,:,2);tDx(i,:,3);tDy(i,:,1);tDy(i,:,2);tDy(i,:,3);tDxx(i,:,1);tDxx(i,:,2);tDxx(i,:,3);....
             tDxy(i,:,1);tDxy(i,:,2);tDxy(i,:,3);tDyy(i,:,1);tDyy(i,:,2);tDyy(i,:,3)];
         
%A(19,:) = nor(1,1)*tDx(i,:,4) + nor(1,2)*tDy(i,:,4);
%A(20,:) = nor(2,1)*tDx(i,:,5) + nor(2,2)*tDy(i,:,5);
%A(21,:) = nor(3,1)*tDx(i,:,6) + nor(3,2)*tDy(i,:,6);

cof = inv(A);
COF((i-1)*base_num+1:i*base_num,:)=cof'; %NT琛屼唬琛ㄨ繖NT涓崟鍏冪i涓熀鍑芥暟鐨勭郴鏁?
end

clear tDx tDy tDxx tDxy tDyy



N = size(node,1);
NT = size(elem,1);

%DG_dof = [(1:NT)', (NT+1:2*NT)', (2*NT+1:3*NT)',(3*NT+1:4*NT)',(4*NT+1:5*NT)',(5*NT+1:6*NT)',(6*NT+1:7*NT)',(7*NT+1:8*NT)',(8*NT+1:9*NT)',(9*NT+1:10*NT)'];
%NDG = max(max(DG_dof));
%BStiff = sparse(3*Ndof,2*NDG);
%DGbase = zeros(DG_base_num,DG_base_num);
%% Information for displacement
DG_base_num = 15;
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
DG_phi(:,10) = lambda2(:,1).^2.*lambda2(:,2).*lambda2(:,3);
DG_phi(:,11) = lambda2(:,1).*lambda2(:,2).^2.*lambda2(:,3);
DG_phi(:,12) = lambda2(:,1).*lambda2(:,2).*lambda2(:,3).^2;
DG_phi(:,13) = lambda2(:,1).^2.*lambda2(:,2).^2;
DG_phi(:,15) = lambda2(:,2).^2.*lambda2(:,3).^2;
DG_phi(:,14) = lambda2(:,3).^2.*lambda2(:,1).^2;
DG_m3 = [3,0,0;0,3,0;0,0,3;2,1,0;2,0,1;0,2,1;1,2,0;1,0,2;0,1,2;2,1,1;1,2,1;1,1,2;2,2,0;2,0,2;0,2,2];
DG_dof = [elem, N+elem, 2*N+elem, (3*N+1:3*N+NT)', (3*N+NT+1:3*N+2*NT)',(3*N+2*NT+1:3*N+3*NT)',(3*N+3*NT+1:3*N+4*NT)',(3*N+4*NT+1:3*N+5*NT)',(3*N+5*NT+1:3*N+6*NT)'];
NDG = max(max(DG_dof));
BStiff = sparse(3*Ndof,2*NDG);

DG_COF = zeros(DG_base_num*NT,DG_base_num);
%DGm = [3,0,0;0,3,0;0,0,3;2,1,0;2,0,1;0,2,1;1,2,0;1,0,2;0,1,2;1,1,1];
DGbase = zeros(DG_base_num,DG_base_num);
DGbase(1,1) = 1;
DGbase(2,2) = 1;
DGbase(3,3) = 1;
for i = 1:DG_base_num
DGbase(10,i) = com_quad(([2,0,0]+DG_m3(i,:))');
DGbase(11,i) = com_quad(([0,2,0]+DG_m3(i,:))');
DGbase(12,i) = com_quad(([0,0,2]+DG_m3(i,:))');
DGbase(13,i) = com_quad(([1,1,0]+DG_m3(i,:))');
DGbase(14,i) = com_quad(([1,0,1]+DG_m3(i,:))');
DGbase(15,i) = com_quad(([0,1,1]+DG_m3(i,:))');
end
for i = 1:NT
DGbase(4,1) = 3*Dlambda(i,1,1);
DGbase(4,4) = Dlambda(i,1,2);
DGbase(4,5) = Dlambda(i,1,3);

DGbase(5,2) = 3*Dlambda(i,1,2); 
DGbase(5,6) = Dlambda(i,1,3);
DGbase(5,7) = Dlambda(i,1,1);

DGbase(6,3) = 3*Dlambda(i,1,3);
DGbase(6,8) = Dlambda(i,1,1);
DGbase(6,9) = Dlambda(i,1,2);

%% y方向的导数
DGbase(7,1) = 3*Dlambda(i,2,1);
DGbase(7,4) = Dlambda(i,2,2);
DGbase(7,5) = Dlambda(i,2,3);

DGbase(8,2) = 3*Dlambda(i,2,2); 
DGbase(8,6) = Dlambda(i,2,3);
DGbase(8,7) = Dlambda(i,2,1);

DGbase(9,3) = 3*Dlambda(i,2,3);
DGbase(9,8) = Dlambda(i,2,1);
DGbase(9,9) = Dlambda(i,2,2);
cof = inv(DGbase);
DG_COF((i-1)*DG_base_num+1:i*DG_base_num,:)=cof'; %NT琛屼唬琛ㄨ繖NT涓崟鍏冪i涓熀鍑芥暟鐨勭郴鏁?
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

Dx(:,10)= 2*lambda2(p,1)*lambda2(p,2)*lambda2(p,3)*Dlambda(:,1,1)+lambda2(p,1)^2*lambda2(p,3)*Dlambda(:,1,2)+....
          lambda2(p,1)^2*lambda2(p,2)*Dlambda(:,1,3);
Dy(:,10)= 2*lambda2(p,1)*lambda2(p,2)*lambda2(p,3)*Dlambda(:,2,1)+lambda2(p,1)^2*lambda2(p,3)*Dlambda(:,2,2)+....
          lambda2(p,1)^2*lambda2(p,2)*Dlambda(:,2,3);

Dx(:,11)= 2*lambda2(p,2)*lambda2(p,1)*lambda2(p,3)*Dlambda(:,1,2)+lambda2(p,2)^2*lambda2(p,3)*Dlambda(:,1,1)+....
          lambda2(p,2)^2*lambda2(p,1)*Dlambda(:,1,3);
Dy(:,11)= 2*lambda2(p,2)*lambda2(p,1)*lambda2(p,3)*Dlambda(:,2,2)+lambda2(p,2)^2*lambda2(p,3)*Dlambda(:,2,1)+....
          lambda2(p,2)^2*lambda2(p,1)*Dlambda(:,2,3);

Dx(:,12)= 2*lambda2(p,3)*lambda2(p,1)*lambda2(p,2)*Dlambda(:,1,3)+lambda2(p,3)^2*lambda2(p,2)*Dlambda(:,1,1)+....
          lambda2(p,3)^2*lambda2(p,1)*Dlambda(:,1,2);
Dy(:,12)= 2*lambda2(p,3)*lambda2(p,1)*lambda2(p,2)*Dlambda(:,2,3)+lambda2(p,3)^2*lambda2(p,2)*Dlambda(:,2,1)+....
          lambda2(p,3)^2*lambda2(p,1)*Dlambda(:,2,2);

Dx(:,13)= 2*lambda2(p,1)*lambda2(p,2)^2*Dlambda(:,1,1)+2*lambda2(p,1)^2*lambda2(p,2)*Dlambda(:,1,2);     
Dy(:,13)= 2*lambda2(p,1)*lambda2(p,2)^2*Dlambda(:,2,1)+2*lambda2(p,1)^2*lambda2(p,2)*Dlambda(:,2,2);         
      
Dx(:,14)= 3*lambda2(p,1)^2*lambda2(p,2)^2*Dlambda(:,1,1)+2*lambda2(p,1)^3*lambda2(p,2)*Dlambda(:,1,2);     
Dy(:,14)= 3*lambda2(p,1)^2*lambda2(p,2)^2*Dlambda(:,2,1)+2*lambda2(p,1)^3*lambda2(p,2)*Dlambda(:,2,2); 

Dx(:,15)= 2*lambda2(p,2)*lambda2(p,3)^2*Dlambda(:,1,2)+2*lambda2(p,2)^2*lambda2(p,3)*Dlambda(:,1,3);     
Dy(:,15)= 2*lambda2(p,2)*lambda2(p,3)^2*Dlambda(:,2,2)+2*lambda2(p,2)^2*lambda2(p,3)*Dlambda(:,2,3);         
      
Dx(:,16)= 3*lambda2(p,2)^2*lambda2(p,3)^2*Dlambda(:,1,2)+2*lambda2(p,2)^3*lambda2(p,3)*Dlambda(:,1,3);     
Dy(:,16)= 3*lambda2(p,2)^2*lambda2(p,3)^2*Dlambda(:,2,2)+2*lambda2(p,2)^3*lambda2(p,3)*Dlambda(:,2,3);

Dx(:,17)= 2*lambda2(p,3)*lambda2(p,1)^2*Dlambda(:,1,3)+2*lambda2(p,3)^2*lambda2(p,1)*Dlambda(:,1,1);     
Dy(:,17)= 2*lambda2(p,3)*lambda2(p,1)^2*Dlambda(:,2,3)+2*lambda2(p,3)^2*lambda2(p,1)*Dlambda(:,2,1);         
      
Dx(:,18)= 3*lambda2(p,3)^2*lambda2(p,1)^2*Dlambda(:,1,3)+2*lambda2(p,3)^3*lambda2(p,1)*Dlambda(:,1,1);     
Dy(:,18)= 3*lambda2(p,3)^2*lambda2(p,1)^2*Dlambda(:,2,3)+2*lambda2(p,3)^3*lambda2(p,1)*Dlambda(:,2,1);
L = [1,2,3];
Dx(:,19) = 2*lambda2(p,L(1))*lambda2(p,L(2))^2*lambda2(p,L(3))*Dlambda(:,1,L(1))+2*lambda2(p,L(1))^2*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,1,L(2))....
             + lambda2(p,L(1))^2*lambda2(p,L(2))^2*Dlambda(:,1,L(3));
Dy(:,19) = 2*lambda2(p,L(1))*lambda2(p,L(2))^2*lambda2(p,L(3))*Dlambda(:,2,L(1))+2*lambda2(p,L(1))^2*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,2,L(2))....
             +  lambda2(p,L(1))^2*lambda2(p,L(2))^2*Dlambda(:,2,L(3)); 
L = [3,1,2];
Dx(:,20) = 2*lambda2(p,L(1))*lambda2(p,L(2))^2*lambda2(p,L(3))*Dlambda(:,1,L(1))+2*lambda2(p,L(1))^2*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,1,L(2))....
             + lambda2(p,L(1))^2*lambda2(p,L(2))^2*Dlambda(:,1,L(3));
Dy(:,20) = 2*lambda2(p,L(1))*lambda2(p,L(2))^2*lambda2(p,L(3))*Dlambda(:,2,L(1))+2*lambda2(p,L(1))^2*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,2,L(2))....
             +  lambda2(p,L(1))^2*lambda2(p,L(2))^2*Dlambda(:,2,L(3));          
L = [2,3,1];
Dx(:,21) = 2*lambda2(p,L(1))*lambda2(p,L(2))^2*lambda2(p,L(3))*Dlambda(:,1,L(1))+2*lambda2(p,L(1))^2*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,1,L(2))....
             + lambda2(p,L(1))^2*lambda2(p,L(2))^2*Dlambda(:,1,L(3));
Dy(:,21) = 2*lambda2(p,L(1))*lambda2(p,L(2))^2*lambda2(p,L(3))*Dlambda(:,2,L(1))+2*lambda2(p,L(1))^2*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,2,L(2))....
             +  lambda2(p,L(1))^2*lambda2(p,L(2))^2*Dlambda(:,2,L(3));         

    
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
%% Fast assembly of BStiff
% The old implementation assembled BStiff inside the quadrature loop by many
% small sparse additions. Here we first accumulate the elementwise local
% integrals, and then assemble the sparse matrix after the quadrature loop.
%
% Blocks:
%   sigma11 -- u1 : d_x Phi_i * Psi_j
%   sigma22 -- u2 : d_y Phi_i * Psi_j
%   sigma12 -- u1 : d_y Phi_i * Psi_j
%   sigma12 -- u2 : d_x Phi_i * Psi_j

Dx = zeros(NT,base_num);
Dy = zeros(NT,base_num);
Bx = zeros(NT,base_num*DG_base_num);  % \int_T d_x Phi_i * Psi_j
By = zeros(NT,base_num*DG_base_num);  % \int_T d_y Phi_i * Psi_j

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

Dx(:,10)= 2*lambda2(p,1)*lambda2(p,2)*lambda2(p,3)*Dlambda(:,1,1)+lambda2(p,1)^2*lambda2(p,3)*Dlambda(:,1,2)+....
          lambda2(p,1)^2*lambda2(p,2)*Dlambda(:,1,3);
Dy(:,10)= 2*lambda2(p,1)*lambda2(p,2)*lambda2(p,3)*Dlambda(:,2,1)+lambda2(p,1)^2*lambda2(p,3)*Dlambda(:,2,2)+....
          lambda2(p,1)^2*lambda2(p,2)*Dlambda(:,2,3);

Dx(:,11)= 2*lambda2(p,2)*lambda2(p,1)*lambda2(p,3)*Dlambda(:,1,2)+lambda2(p,2)^2*lambda2(p,3)*Dlambda(:,1,1)+....
          lambda2(p,2)^2*lambda2(p,1)*Dlambda(:,1,3);
Dy(:,11)= 2*lambda2(p,2)*lambda2(p,1)*lambda2(p,3)*Dlambda(:,2,2)+lambda2(p,2)^2*lambda2(p,3)*Dlambda(:,2,1)+....
          lambda2(p,2)^2*lambda2(p,1)*Dlambda(:,2,3);

Dx(:,12)= 2*lambda2(p,3)*lambda2(p,1)*lambda2(p,2)*Dlambda(:,1,3)+lambda2(p,3)^2*lambda2(p,2)*Dlambda(:,1,1)+....
          lambda2(p,3)^2*lambda2(p,1)*Dlambda(:,1,2);
Dy(:,12)= 2*lambda2(p,3)*lambda2(p,1)*lambda2(p,2)*Dlambda(:,2,3)+lambda2(p,3)^2*lambda2(p,2)*Dlambda(:,2,1)+....
          lambda2(p,3)^2*lambda2(p,1)*Dlambda(:,2,2);

Dx(:,13)= 2*lambda2(p,1)*lambda2(p,2)^2*Dlambda(:,1,1)+2*lambda2(p,1)^2*lambda2(p,2)*Dlambda(:,1,2);     
Dy(:,13)= 2*lambda2(p,1)*lambda2(p,2)^2*Dlambda(:,2,1)+2*lambda2(p,1)^2*lambda2(p,2)*Dlambda(:,2,2);         
      
Dx(:,14)= 3*lambda2(p,1)^2*lambda2(p,2)^2*Dlambda(:,1,1)+2*lambda2(p,1)^3*lambda2(p,2)*Dlambda(:,1,2);     
Dy(:,14)= 3*lambda2(p,1)^2*lambda2(p,2)^2*Dlambda(:,2,1)+2*lambda2(p,1)^3*lambda2(p,2)*Dlambda(:,2,2); 

Dx(:,15)= 2*lambda2(p,2)*lambda2(p,3)^2*Dlambda(:,1,2)+2*lambda2(p,2)^2*lambda2(p,3)*Dlambda(:,1,3);     
Dy(:,15)= 2*lambda2(p,2)*lambda2(p,3)^2*Dlambda(:,2,2)+2*lambda2(p,2)^2*lambda2(p,3)*Dlambda(:,2,3);         
      
Dx(:,16)= 3*lambda2(p,2)^2*lambda2(p,3)^2*Dlambda(:,1,2)+2*lambda2(p,2)^3*lambda2(p,3)*Dlambda(:,1,3);     
Dy(:,16)= 3*lambda2(p,2)^2*lambda2(p,3)^2*Dlambda(:,2,2)+2*lambda2(p,2)^3*lambda2(p,3)*Dlambda(:,2,3);

Dx(:,17)= 2*lambda2(p,3)*lambda2(p,1)^2*Dlambda(:,1,3)+2*lambda2(p,3)^2*lambda2(p,1)*Dlambda(:,1,1);     
Dy(:,17)= 2*lambda2(p,3)*lambda2(p,1)^2*Dlambda(:,2,3)+2*lambda2(p,3)^2*lambda2(p,1)*Dlambda(:,2,1);         
      
Dx(:,18)= 3*lambda2(p,3)^2*lambda2(p,1)^2*Dlambda(:,1,3)+2*lambda2(p,3)^3*lambda2(p,1)*Dlambda(:,1,1);     
Dy(:,18)= 3*lambda2(p,3)^2*lambda2(p,1)^2*Dlambda(:,2,3)+2*lambda2(p,3)^3*lambda2(p,1)*Dlambda(:,2,1);
L = [1,2,3];
Dx(:,19) = 2*lambda2(p,L(1))*lambda2(p,L(2))^2*lambda2(p,L(3))*Dlambda(:,1,L(1))+2*lambda2(p,L(1))^2*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,1,L(2))....
             + lambda2(p,L(1))^2*lambda2(p,L(2))^2*Dlambda(:,1,L(3));
Dy(:,19) = 2*lambda2(p,L(1))*lambda2(p,L(2))^2*lambda2(p,L(3))*Dlambda(:,2,L(1))+2*lambda2(p,L(1))^2*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,2,L(2))....
             +  lambda2(p,L(1))^2*lambda2(p,L(2))^2*Dlambda(:,2,L(3)); 
L = [3,1,2];
Dx(:,20) = 2*lambda2(p,L(1))*lambda2(p,L(2))^2*lambda2(p,L(3))*Dlambda(:,1,L(1))+2*lambda2(p,L(1))^2*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,1,L(2))....
             + lambda2(p,L(1))^2*lambda2(p,L(2))^2*Dlambda(:,1,L(3));
Dy(:,20) = 2*lambda2(p,L(1))*lambda2(p,L(2))^2*lambda2(p,L(3))*Dlambda(:,2,L(1))+2*lambda2(p,L(1))^2*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,2,L(2))....
             +  lambda2(p,L(1))^2*lambda2(p,L(2))^2*Dlambda(:,2,L(3));          
L = [2,3,1];
Dx(:,21) = 2*lambda2(p,L(1))*lambda2(p,L(2))^2*lambda2(p,L(3))*Dlambda(:,1,L(1))+2*lambda2(p,L(1))^2*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,1,L(2))....
             + lambda2(p,L(1))^2*lambda2(p,L(2))^2*Dlambda(:,1,L(3));
Dy(:,21) = 2*lambda2(p,L(1))*lambda2(p,L(2))^2*lambda2(p,L(3))*Dlambda(:,2,L(1))+2*lambda2(p,L(1))^2*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,2,L(2))....
             +  lambda2(p,L(1))^2*lambda2(p,L(2))^2*Dlambda(:,2,L(3));         

    % Evaluate all stress derivative basis functions at this quadrature point
    Sx = zeros(NT,base_num);
    Sy = zeros(NT,base_num);
    for ii = 1:base_num
        COFii = COF(ii:base_num:end,:);
        Sx(:,ii) = sum(COFii.*Dx,2);
        Sy(:,ii) = sum(COFii.*Dy,2);
    end

    % Evaluate all displacement basis functions at this quadrature point
    Uh = zeros(NT,DG_base_num);
    for jj = 1:DG_base_num
        Uh(:,jj) = DG_COF(jj:DG_base_num:end,:)*DG_phi(p,:)';
    end

    % Accumulate local blocks before multiplying by area
    for ii = 1:base_num
        cols = (ii-1)*DG_base_num + (1:DG_base_num);
        Bx(:,cols) = Bx(:,cols) + weight2(p)*bsxfun(@times,Sx(:,ii),Uh);
        By(:,cols) = By(:,cols) + weight2(p)*bsxfun(@times,Sy(:,ii),Uh);
    end

    clear Sx Sy Uh COFii cols
end

% Multiply the elementwise integrals by the physical element areas
Bx = bsxfun(@times,Bx,area);
By = bsxfun(@times,By,area);

% Assemble BStiff after quadrature. This avoids repeated sparse additions
% inside the quadrature loop.
BStiff = sparse(3*Ndof,2*NDG);
for ii = 1:base_num
    row11 = double(elem2dof(:,ii));
    row22 = Ndof + double(elem2dof(:,ii));
    row12 = 2*Ndof + double(elem2dof(:,ii));

    for jj = 1:DG_base_num
        idx = (ii-1)*DG_base_num + jj;
        col1 = double(DG_dof(:,jj));
        col2 = NDG + double(DG_dof(:,jj));

        BStiff = BStiff + sparse(row11,col1,Bx(:,idx),3*Ndof,2*NDG);
        BStiff = BStiff + sparse(row22,col2,By(:,idx),3*Ndof,2*NDG);
        BStiff = BStiff + sparse(row12,col1,By(:,idx),3*Ndof,2*NDG);
        BStiff = BStiff + sparse(row12,col2,Bx(:,idx),3*Ndof,2*NDG);
    end
end

% Clear large temporary arrays for memory saving
clear Bx By Dx Dy row11 row22 row12 col1 col2 idx ii jj p
%Stiff = Stiff + (1/hsize^5)*pA;
%% 鍙崇椤?        
%big_Stiff = [iota*Stiff+M, sparse(Ndof,Ndof),sparse(Ndof,Ndof);sparse(Ndof,Ndof),iota*Stiff+M,sparse(Ndof,Ndof);sparse(Ndof,Ndof),sparse(Ndof,Ndof),2*iota*Stiff+2*M];%blkdiag(Stiff,Stiff);
beta_m = (2*mu+La)/(4*mu*(La+mu));
beta_G = -La/(4*mu*(La+mu));
big_Stiff = [beta_m*(iota^2*Stiff+M),beta_G*(iota^2*Stiff+M) ,sparse(Ndof,Ndof);beta_G*(iota^2*Stiff+M),beta_m*(iota^2*Stiff+M),sparse(Ndof,Ndof);sparse(Ndof,Ndof),sparse(Ndof,Ndof),(iota^2*Stiff+M)/mu];%blkdiag(Stiff,Stiff);
All_Stiff = [big_Stiff,BStiff;BStiff',sparse(2*NDG,2*NDG)];
clear Stiff M BStiff Bij11 Bij22 Bij31 Bij32 
%% Right hand side
Bt = zeros(NT,DG_base_num);
Bt1 = Bt;
 for p  =1:Quad2        
  pxy = lambda2(p,1)*node(elem(:,1),:)+lambda2(p,2)*node(elem(:,2),:)+lambda2(p,3)*node(elem(:,3),:);   
 % [f1,f2]=Ipf_Straingradient(pxy,mu,La,ell);
 [f,f1] = layer_Ipf_LSG(pxy,La,mu);
 %[f,f1] = La_Ipf_LSG(pxy,La,mu);

 %f = ones(NT,1);
%f1 = ones(NT,1);
% for i = 1:base_num
%   Bt(:,i)=Bt(:,i)+weight2(p)*(COF(i:base_num:end,:)*(phi(p,:)')).*f;  
%   Bt1(:,i)=Bt1(:,i)+weight2(p)*(COF(i:base_num:end,:)*(phi(p,:)')).*f1; 
% end
  
for i = 1:DG_base_num
  Bt(:,i)=Bt(:,i)+weight2(p)*(DG_COF(i:DG_base_num:end,:)*DG_phi(p,:)').*f;  
  Bt1(:,i)=Bt1(:,i)+weight2(p)*(DG_COF(i:DG_base_num:end,:)*DG_phi(p,:)').*f1; 
end

 end
    Bt = Bt.*repmat(area,1,DG_base_num);
    Bt1 = Bt1.*repmat(area,1,DG_base_num);
    b = accumarray(DG_dof(:),Bt(:),[NDG 1]);
    b1 = accumarray(DG_dof(:),Bt1(:),[NDG 1]);
    right = [zeros(3*Ndof,1);b;b1];
  %  right = [zeros(3*Ndof,1);ones(2*NDG,1)];
    solut = All_Stiff\right;

solut1 = solut(3*Ndof+1:3*Ndof+NDG);
solut2 = solut(3*Ndof+1+NDG:3*Ndof+2*NDG);
solut_s11 = solut(1:Ndof);
solut_s22 = solut(1+Ndof:2*Ndof);
solut_s12 = solut(2*Ndof+1:3*Ndof);
intu1 = solut1(DG_dof);
intu2 = solut2(DG_dof);
int_s11 = solut_s11(elem2dof);
int_s22 = solut_s22(elem2dof);
int_s12 = solut_s12(elem2dof);


%intu = uh(elem2dof);
error = zeros(NT,1);
enor = zeros(NT,1);
error_s = zeros(NT,1);
for p = 1:Quad2
    ipinteru1 = zeros(NT,DG_base_num);
    ipinteru2 = ipinteru1;
    ipinters11 = zeros(NT,base_num);
    ipinters22 = ipinters11;
    ipinters12 = ipinters11;
  for i = 1:DG_base_num
   ipinteru1(:,i)= ipinteru1(:,i)+intu1(:,i).*(DG_COF(i:DG_base_num:end,:)*DG_phi(p,:)');
   ipinteru2(:,i)= ipinteru2(:,i)+intu2(:,i).*(DG_COF(i:DG_base_num:end,:)*DG_phi(p,:)');
  end  
  for i = 1:base_num
    ipinters11(:,i)= ipinters11(:,i)+int_s11(:,i).*(COF(i:base_num:end,:)*(phi(p,:)'));
    ipinters12(:,i)= ipinters12(:,i)+int_s12(:,i).*(COF(i:base_num:end,:)*(phi(p,:)'));
    ipinters22(:,i)= ipinters22(:,i)+int_s22(:,i).*(COF(i:base_num:end,:)*(phi(p,:)'));
  end
  pxy = lambda2(p,1)*node(elem(:,1),:)+lambda2(p,2)*node(elem(:,2),:)+lambda2(p,3)*node(elem(:,3),:);  
  [ru1,ru2] = layer_Ipr_LSG(pxy, La);
  [sigma11,sigma22,sigma12] = layer_sigma_Ipr_LSG(pxy,La,mu);
  ipiu1 = sum(ipinteru1,2);
  ipiu2 = sum(ipinteru2,2);
  ipi_s11 = sum(ipinters11,2);
  ipi_s22 = sum(ipinters22,2);
  ipi_s12 = sum(ipinters12,2);
%  enor = enor + weight2(p)*ipi.^2;
 % error = error + weight2(p)*((ru1-ipiu1).^2 + (ru2-ipiu2).^2  );%+weight2(p)*(ru2-ipi1).^2;
  %error_s = error_s + weight2(p)*((sigma11-ipi_s11).^2 + (sigma12-ipi_s12).^2  + (sigma22-ipi_s22).^2   );
    error = error + weight2(p)*((ru1-ipiu1).^2 + (ru2-ipiu2).^2);%+weight2(p)*(ru2-ipi1).^2;
  error_s = error_s + weight2(p)*((sigma11-ipi_s11).^2 + (sigma22-ipi_s22).^2  + 2*(sigma12-ipi_s12).^2 );
%   error = error + weight2(p)*(ipi).^2;
%   error_s = error_s + weight2(p)*(ipi_s).^2;  
end

error2 = zeros(NT,1);
Dx = zeros(NT,base_num);
Dy = zeros(NT,base_num);


    Dfx11 = zeros(NT,base_num);
    Dfy11 = zeros(NT,base_num);
    Dfx12 = zeros(NT,base_num);
    Dfy12 = zeros(NT,base_num);  
    Dfx22 = zeros(NT,base_num);
    Dfy22 = zeros(NT,base_num);  
    
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

Dx(:,10)= 2*lambda2(p,1)*lambda2(p,2)*lambda2(p,3)*Dlambda(:,1,1)+lambda2(p,1)^2*lambda2(p,3)*Dlambda(:,1,2)+....
          lambda2(p,1)^2*lambda2(p,2)*Dlambda(:,1,3);
Dy(:,10)= 2*lambda2(p,1)*lambda2(p,2)*lambda2(p,3)*Dlambda(:,2,1)+lambda2(p,1)^2*lambda2(p,3)*Dlambda(:,2,2)+....
          lambda2(p,1)^2*lambda2(p,2)*Dlambda(:,2,3);

Dx(:,11)= 2*lambda2(p,2)*lambda2(p,1)*lambda2(p,3)*Dlambda(:,1,2)+lambda2(p,2)^2*lambda2(p,3)*Dlambda(:,1,1)+....
          lambda2(p,2)^2*lambda2(p,1)*Dlambda(:,1,3);
Dy(:,11)= 2*lambda2(p,2)*lambda2(p,1)*lambda2(p,3)*Dlambda(:,2,2)+lambda2(p,2)^2*lambda2(p,3)*Dlambda(:,2,1)+....
          lambda2(p,2)^2*lambda2(p,1)*Dlambda(:,2,3);

Dx(:,12)= 2*lambda2(p,3)*lambda2(p,1)*lambda2(p,2)*Dlambda(:,1,3)+lambda2(p,3)^2*lambda2(p,2)*Dlambda(:,1,1)+....
          lambda2(p,3)^2*lambda2(p,1)*Dlambda(:,1,2);
Dy(:,12)= 2*lambda2(p,3)*lambda2(p,1)*lambda2(p,2)*Dlambda(:,2,3)+lambda2(p,3)^2*lambda2(p,2)*Dlambda(:,2,1)+....
          lambda2(p,3)^2*lambda2(p,1)*Dlambda(:,2,2);

Dx(:,13)= 2*lambda2(p,1)*lambda2(p,2)^2*Dlambda(:,1,1)+2*lambda2(p,1)^2*lambda2(p,2)*Dlambda(:,1,2);     
Dy(:,13)= 2*lambda2(p,1)*lambda2(p,2)^2*Dlambda(:,2,1)+2*lambda2(p,1)^2*lambda2(p,2)*Dlambda(:,2,2);         
      
Dx(:,14)= 3*lambda2(p,1)^2*lambda2(p,2)^2*Dlambda(:,1,1)+2*lambda2(p,1)^3*lambda2(p,2)*Dlambda(:,1,2);     
Dy(:,14)= 3*lambda2(p,1)^2*lambda2(p,2)^2*Dlambda(:,2,1)+2*lambda2(p,1)^3*lambda2(p,2)*Dlambda(:,2,2); 

Dx(:,15)= 2*lambda2(p,2)*lambda2(p,3)^2*Dlambda(:,1,2)+2*lambda2(p,2)^2*lambda2(p,3)*Dlambda(:,1,3);     
Dy(:,15)= 2*lambda2(p,2)*lambda2(p,3)^2*Dlambda(:,2,2)+2*lambda2(p,2)^2*lambda2(p,3)*Dlambda(:,2,3);         
      
Dx(:,16)= 3*lambda2(p,2)^2*lambda2(p,3)^2*Dlambda(:,1,2)+2*lambda2(p,2)^3*lambda2(p,3)*Dlambda(:,1,3);     
Dy(:,16)= 3*lambda2(p,2)^2*lambda2(p,3)^2*Dlambda(:,2,2)+2*lambda2(p,2)^3*lambda2(p,3)*Dlambda(:,2,3);

Dx(:,17)= 2*lambda2(p,3)*lambda2(p,1)^2*Dlambda(:,1,3)+2*lambda2(p,3)^2*lambda2(p,1)*Dlambda(:,1,1);     
Dy(:,17)= 2*lambda2(p,3)*lambda2(p,1)^2*Dlambda(:,2,3)+2*lambda2(p,3)^2*lambda2(p,1)*Dlambda(:,2,1);         
      
Dx(:,18)= 3*lambda2(p,3)^2*lambda2(p,1)^2*Dlambda(:,1,3)+2*lambda2(p,3)^3*lambda2(p,1)*Dlambda(:,1,1);     
Dy(:,18)= 3*lambda2(p,3)^2*lambda2(p,1)^2*Dlambda(:,2,3)+2*lambda2(p,3)^3*lambda2(p,1)*Dlambda(:,2,1);
L = [1,2,3];
Dx(:,19) = 2*lambda2(p,L(1))*lambda2(p,L(2))^2*lambda2(p,L(3))*Dlambda(:,1,L(1))+2*lambda2(p,L(1))^2*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,1,L(2))....
             + lambda2(p,L(1))^2*lambda2(p,L(2))^2*Dlambda(:,1,L(3));
Dy(:,19) = 2*lambda2(p,L(1))*lambda2(p,L(2))^2*lambda2(p,L(3))*Dlambda(:,2,L(1))+2*lambda2(p,L(1))^2*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,2,L(2))....
             +  lambda2(p,L(1))^2*lambda2(p,L(2))^2*Dlambda(:,2,L(3)); 
L = [3,1,2];
Dx(:,20) = 2*lambda2(p,L(1))*lambda2(p,L(2))^2*lambda2(p,L(3))*Dlambda(:,1,L(1))+2*lambda2(p,L(1))^2*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,1,L(2))....
             + lambda2(p,L(1))^2*lambda2(p,L(2))^2*Dlambda(:,1,L(3));
Dy(:,20) = 2*lambda2(p,L(1))*lambda2(p,L(2))^2*lambda2(p,L(3))*Dlambda(:,2,L(1))+2*lambda2(p,L(1))^2*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,2,L(2))....
             +  lambda2(p,L(1))^2*lambda2(p,L(2))^2*Dlambda(:,2,L(3));          
L = [2,3,1];
Dx(:,21) = 2*lambda2(p,L(1))*lambda2(p,L(2))^2*lambda2(p,L(3))*Dlambda(:,1,L(1))+2*lambda2(p,L(1))^2*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,1,L(2))....
             + lambda2(p,L(1))^2*lambda2(p,L(2))^2*Dlambda(:,1,L(3));
Dy(:,21) = 2*lambda2(p,L(1))*lambda2(p,L(2))^2*lambda2(p,L(3))*Dlambda(:,2,L(1))+2*lambda2(p,L(1))^2*lambda2(p,L(2))*lambda2(p,L(3))*Dlambda(:,2,L(2))....
             +  lambda2(p,L(1))^2*lambda2(p,L(2))^2*Dlambda(:,2,L(3));    
for i = 1:base_num
    Dfx11(:,i) =  int_s11(:,i).*dot(COF(i:base_num:end,:),Dx,2);
    Dfy11(:,i) =  int_s11(:,i).*dot(COF(i:base_num:end,:),Dy,2);
    Dfx12(:,i) =  int_s12(:,i).*dot(COF(i:base_num:end,:),Dx,2);
    Dfy12(:,i) =  int_s12(:,i).*dot(COF(i:base_num:end,:),Dy,2);
    Dfx22(:,i) =  int_s22(:,i).*dot(COF(i:base_num:end,:),Dx,2);
    Dfy22(:,i) =  int_s22(:,i).*dot(COF(i:base_num:end,:),Dy,2);

end
  pxy = lambda2(p,1)*node(elem(:,1),:)+lambda2(p,2)*node(elem(:,2),:)+lambda2(p,3)*node(elem(:,3),:);  
  iDfx11 = sum(Dfx11,2);
  iDfy11 = sum(Dfy11,2);
  iDfx12 = sum(Dfx12,2);
  iDfy12 = sum(Dfy12,2);
  iDfx22 = sum(Dfx22,2);
  iDfy22 = sum(Dfy22,2);
 % iDfxy = sum(Dfxy,2);
 [ s11y, s22y , s12y ] = layer_sigmay_Ipr_LSG( pxy ,La, mu );
 [ s11x, s22x , s12x ] = layer_sigmax_Ipr_LSG( pxy ,La, mu );
  
  
error2 = error2 + weight2(p)*((iDfx11+iDfy12-s11x-s12y).^2 +(iDfx12+iDfy22-s12x-s22y).^2  );
error2 = error2 +  iota^2*weight2(p) * ( ...
      (iDfx11 - s11x).^2 + (iDfy11 - s11y).^2 ...
    + (iDfx22 - s22x).^2 + (iDfy22 - s22y).^2 ...
    + 2*((iDfx12 - s12x).^2 + (iDfy12 - s12y).^2) );
%error2 = error2 + ell^2*weight2(p)*((iDfxx-fxx).^2+(iDfxy-fxy).^2+(iDfyy-fyy).^2 +(iDfxx1-fxx1).^2+(iDfxy1-fxy1).^2+(iDfyy1-fyy1).^2);
%error2 = error2 + weight2(p)*(iDfx-fx).^2;
 end    



error = error.*area;
error2 = error2.*area;
enor = enor.*area;
enor = sqrt(sum(enor));
error1 = sum(error);
error1 = sqrt(error1);

error_s=error_s.*area;
errors = sum(error_s);
errors = sqrt(errors);
error_e = sqrt(sum(error2));
error_infty = max(abs(solut(1+5*size(node,1):6*size(node,1))));
diverror = sqrt(sum(error2))+errors;
disp(['Hsize:=' num2str(hsize)]);
%disp(['L2 Error of u1:=' num2str(error1)]);
%disp(['L2 Error of sigma_11:=' num2str(errors)]);
%disp(['div Error of sigma_11:=' num2str(diverror)]);
%disp(['infty Error of sigma_11:=' num2str(error_infty)]);
fprintf('iota = %.1e, energy norm error = %.16e\n', iota, diverror+error1);
% if I_k==0
% save('solutiota1.mat','solut');
% else if I_k == 2
%       save('solutiota10e2.mat','solut');  
%     else if I_k == 4 
%       save('solutiota10e4.mat','solut');
%         else if I_k == 6
%       save('solutiota10e6.mat','solut');
%             else if I_k == 8
%       save('solutiota10e8.mat','solut');  
%                 end
%             end
%         end
%     end
% end
% clear

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


