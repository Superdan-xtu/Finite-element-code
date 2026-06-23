%[node,elem] = squaremesh([0,1,0,1],1/128);
hsize = 1/64;
[X,Y] = meshgrid(0:1/64:1, 0:1/64:1);
% refinesolut1 = load('solutiota1new.mat');
% solut = refinesolut1.solut;
% refinemesh1 = load('elemla100000iota1.mat');
% Relem = refinemesh1.elem;
% refinenode1 = load('nodela100000iota1.mat');
% Rnode = refinenode1.node;
% refinestresscof = load('COFla100000iota1.mat');
% COF = refinestresscof.COF;
% refineDGcof = load('DGCOFla100000iota1.mat');
% DG_COF = refineDGcof.DG_COF;
% NT = size(elem,1);
% N = size(node,1);
% [elem2dof1,elem2edge,edge,bdDof,freeDof] = C2dof(elem,node); 
% Ndof = double(max(elem2dof1(:)));
% elem2dof = [elem2dof1];
% DG_dof = [elem, N+elem, 2*N+elem, (3*N+1:3*N+NT)', (3*N+NT+1:3*N+2*NT)',(3*N+2*NT+1:3*N+3*NT)',(3*N+3*NT+1:3*N+4*NT)',(3*N+4*NT+1:3*N+5*NT)',(3*N+5*NT+1:3*N+6*NT)'];
% NDG = max(max(DG_dof));
%  [node,elem] = squaremesh([0,1,0,1],hsize); 
 sigma11temp = solut(1+5*size(node,1):6*size(node,1));
% %sigma11temp = solut(3*Ndof+size(node,1)+1:3*Ndof+2*size(node,1));
sigma11= reshape(sigma11temp,65,65);
figure
surf(X,Y,sigma11);
% fL2 = zeros(NT,1);
% for p = 1:Quad2
%    pxy = lambda2(p,1)*node(elem(:,1),:)+lambda2(p,2)*node(elem(:,2),:)+lambda2(p,3)*node(elem(:,3),:);     
%      [f,f1] = layer_Ipf_LSG(pxy,La,mu);
%      fL2 = fL2 + weight2(p)*(f.^2);
% end
%normL2f = sqrt(sum(fL2.*area)); %  3.404901530440675
%%%%%%%%%
% order2 = 9;
% [lambda2,weight2] = quadpts(order2);
% Quad2 = size(lambda2,1);
% solut1 = solut(3*Ndof+1:3*Ndof+NDG);
% solut2 = solut(3*Ndof+1+NDG:3*Ndof+2*NDG);
% solut_s11 = solut(1:Ndof);
% solut_s22 = solut(1+Ndof:2*Ndof);
% solut_s12 = solut(2*Ndof+1:3*Ndof);
% intu1 = solut1(DG_dof);
% intu2 = solut2(DG_dof);
% int_s11 = solut_s11(elem2dof);
% int_s22 = solut_s22(elem2dof);
% int_s12 = solut_s12(elem2dof);
% 
% 
% %intu = uh(elem2dof);
% error = zeros(NT,1);
% enor = zeros(NT,1);
% error_s = zeros(NT,1);
% for p = 1:Quad2
% %    ipinteru1 = zeros(NT,DG_base_num);
% %    ipinteru2 = ipinteru1;
%     ipinters11 = zeros(NT,base_num);
%     ipinters22 = ipinters11;
%     ipinters12 = ipinters11;
% %   for i = 1:DG_base_num
% %    ipinteru1(:,i)= ipinteru1(:,i)+intu1(:,i).*(DG_COF(i:DG_base_num:end,:)*DG_phi(p,:)');
% %    ipinteru2(:,i)= ipinteru2(:,i)+intu2(:,i).*(DG_COF(i:DG_base_num:end,:)*DG_phi(p,:)');
% %   end  
%   for i = 1:base_num
%     ipinters11(:,i)= ipinters11(:,i)+int_s11(:,i).*(COF(i:base_num:end,:)*(phi(p,:)'));
%     ipinters12(:,i)= ipinters12(:,i)+int_s12(:,i).*(COF(i:base_num:end,:)*(phi(p,:)'));
%     ipinters22(:,i)= ipinters22(:,i)+int_s22(:,i).*(COF(i:base_num:end,:)*(phi(p,:)'));
%   end
%   pxy = lambda2(p,1)*node(elem(:,1),:)+lambda2(p,2)*node(elem(:,2),:)+lambda2(p,3)*node(elem(:,3),:);  
%   [ru1,ru2] = layer_Ipr_LSG(pxy, La);
%   %[ f1, f2, f3 ] = La_Ipr_sigma_LSG( X,La,mu )
%   [sigma11,sigma22,sigma12] = La_Ipr_sigma_LSG(pxy,La,mu);
% %  ipiu1 = sum(ipinteru1,2);
% %  ipiu2 = sum(ipinteru2,2);
%   ipi_s11 = sum(ipinters11,2);
%   ipi_s22 = sum(ipinters22,2);
%   ipi_s12 = sum(ipinters12,2);
% %  enor = enor + weight2(p)*ipi.^2;
%  % error = error + weight2(p)*((ru1-ipiu1).^2 + (ru2-ipiu2).^2  );%+weight2(p)*(ru2-ipi1).^2;
%   %error_s = error_s + weight2(p)*((sigma11-ipi_s11).^2 + (sigma12-ipi_s12).^2  + (sigma22-ipi_s22).^2   );
%  %   error = error + weight2(p)*((ru1-ipiu1).^2 );%+weight2(p)*(ru2-ipi1).^2;
%   error_s = error_s + weight2(p)*((sigma11-ipi_s11).^2  );
% %   error = error + weight2(p)*(ipi).^2;
% %   error_s = error_s + weight2(p)*(ipi_s).^2;  
% end