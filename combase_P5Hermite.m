function [ A] = combase_P5Hermite( Dlambda)
%IP_BASE Summary of this function goes here
%   Detailed explanation goes here
k = 21;
A = zeros(k,k);
%% 涓変釜椤剁偣鐨勮嚜鐢卞害
A(1,1)=1;
%A(1,10)=1/10;
A(2,2)=1;
%A(2,10)=1/10;
A(3,3)=1;
%A(3,10)=1/10;

%% x方向的导数
A(4,1) = 5*Dlambda(1,1);
A(4,4) = Dlambda(1,2);
A(4,5) = Dlambda(1,3);

A(5,2) = 5*Dlambda(1,2); 
A(5,6) = Dlambda(1,3);
A(5,7) = Dlambda(1,1);

A(6,3) = 5*Dlambda(1,3);
A(6,8) = Dlambda(1,1);
A(6,9) = Dlambda(1,2);

%% y方向的导数
A(7,1) = 5*Dlambda(2,1);
A(7,4) = Dlambda(2,2);
A(7,5) = Dlambda(2,3);

A(8,2) = 5*Dlambda(2,2); 
A(8,6) = Dlambda(2,3);
A(8,7) = Dlambda(2,1);

A(9,3) = 5*Dlambda(2,3);
A(9,8) = Dlambda(2,1);
A(9,9) = Dlambda(2,2);


%% 第一条边
%   A(10,1)=0; 
%   
%   A(10,2)=3*QA(1,1)*(Dlambda(:,2)'*nor(:,1));
%   
%   A(10,3)=3*QA(2,1)*(Dlambda(:,3)'*nor(:,1));
%   %A(5,3)=3*QA(2,2)*(Dlambda(:,3)'*nor(:,1));
%   
%   A(10,4)=0;
%   
%   A(10,5)=0;
%  
%   A(10,6)=2*QA(3,1)*(Dlambda(:,2)'*nor(:,1))+QA(1,1)*(Dlambda(:,3)'*nor(:,1));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
%   %A(5,6)=2*QA(3,2)*(Dlambda(:,2)'*nor(:,1))+QA(1,2)*(Dlambda(:,3)'*nor(:,1));
%   
%   A(10,7)=QA(1,1)*(Dlambda(:,1)'*nor(:,1));
%   %A(5,7)=QA(1,2)*(Dlambda(:,1)'*nor(:,1));
%   
%   A(10,8)=QA(2,1)*(Dlambda(:,1)'*nor(:,1));
%   %A(5,8)=QA(2,2)*(Dlambda(:,1)'*nor(:,1));
%   
%   A(10,9)=QA(2,1)*(Dlambda(:,2)'*nor(:,1))+2*QA(3,1)*(Dlambda(:,3)'*nor(:,1));
%   %A(5,9)=QA(2,2)*(Dlambda(:,2)'*nor(:,1))+2*QA(3,2)*(Dlambda(:,3)'*nor(:,1));
%   
%   %A(10,10)= -Dlambda(:,1)'*nor(:,1); %%% 前面系数是负的？
%   %A(5,10)=0;
%   A(10,11) = QA(4,1)*Dlambda(:,1)'*nor(:,1);
%   A(10,12) = QA(5,1)*Dlambda(:,1)'*nor(:,1);
%   
%   A(10,15) = 2*QA(5,1)*Dlambda(:,2)'*nor(:,1) + 2*QA(4,1)*Dlambda(:,3)'*nor(:,1);
% 
%   
%   
%   %% 第二条边
%   A(11,2)=0; 
%   %A(7,2)=0;
%   
%   A(11,1)=3*QA(2,1)*(Dlambda(:,1)'*nor(:,2));
%   %A(7,1)=3*QA(2,2)*(Dlambda(:,1)'*nor(:,2));
%   
%   A(11,3)=3*QA(1,1)*(Dlambda(:,3)'*nor(:,2));
%   %A(7,3)=3*QA(1,2)*(Dlambda(:,3)'*nor(:,2));
%   
%   A(11,6)=0;
%   %A(7,6)=0;
%   
%   A(11,7)=0;
%   %A(7,7)=0;
%   
%   A(11,8)=2*QA(3,1)*(Dlambda(:,3)'*nor(:,2))+QA(1,1)*(Dlambda(:,1)'*nor(:,2));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
%   %A(7,8)=2*QA(3,2)*(Dlambda(:,3)'*nor(:,2))+QA(1,2)*(Dlambda(:,1)'*nor(:,2));
%   
%   A(11,9)=QA(1,1)*(Dlambda(:,2)'*nor(:,2));
%   %A(7,9)=QA(1,2)*(Dlambda(:,2)'*nor(:,2));
%   
%   A(11,4)=QA(2,1)*(Dlambda(:,2)'*nor(:,2));
%   %A(7,4)=QA(2,2)*(Dlambda(:,2)'*nor(:,2));
%   
%   A(11,5)=QA(2,1)*(Dlambda(:,3)'*nor(:,2))+2*QA(3,1)*(Dlambda(:,1)'*nor(:,2));
%   %A(7,5)=QA(2,2)*(Dlambda(:,3)'*nor(:,2))+2*QA(3,2)*(Dlambda(:,1)'*nor(:,2));
%   
%   %A(11,10)=QA(5,1)*(Dlambda(:,2)'*nor(:,2));
%   %A(7,10)=QA(5,2)*(Dlambda(:,2)'*nor(:,2));
%   
%   %A(11,11)=-Dlambda(:,2)'*nor(:,2);
% %  A(7,11)=0;
%  A(11,10) = QA(5,1)*Dlambda(:,2)'*nor(:,2);
%  A(11,12) = QA(4,1)*Dlambda(:,2)'*nor(:,2);
%  A(11,14) = 2*QA(4,1)*Dlambda(:,1)'*nor(:,2) + 2*QA(5,1)*Dlambda(:,3)'*nor(:,2);
%  %% 第三条边
%   
%   A(12,3)=0; 
%   %A(9,3)=0;
%   
%   A(12,1)=3*QA(1,1)*(Dlambda(:,1)'*nor(:,3));
%  % A(9,1)=3*QA(1,2)*(Dlambda(:,1)'*nor(:,3));
%   
%   A(12,2)=3*QA(2,1)*(Dlambda(:,2)'*nor(:,3));
%   %A(9,2)=3*QA(2,2)*(Dlambda(:,2)'*nor(:,3));
%   
%   A(12,8)=0;
%   %A(9,8)=0;
%   
%   A(12,9)=0;
%   %A(9,9)=0;
%   
%   A(12,4)=2*QA(3,1)*(Dlambda(:,1)'*nor(:,3))+QA(1,1)*(Dlambda(:,2)'*nor(:,3));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
%   %A(9,4)=2*QA(3,2)*(Dlambda(:,1)'*nor(:,3))+QA(1,2)*(Dlambda(:,2)'*nor(:,3));
%   
%   A(12,5)=QA(1,1)*(Dlambda(:,3)'*nor(:,3));
%   %A(9,5)=QA(1,2)*(Dlambda(:,3)'*nor(:,3));
%   
%   A(12,6)=QA(2,1)*(Dlambda(:,3)'*nor(:,3));
%   %A(9,6)=QA(2,2)*(Dlambda(:,3)'*nor(:,3));
%   
%   A(12,7)=QA(2,1)*(Dlambda(:,1)'*nor(:,3))+2*QA(3,1)*(Dlambda(:,2)'*nor(:,3));
%   %A(9,7)=QA(2,2)*(Dlambda(:,1)'*nor(:,3))+2*QA(3,2)*(Dlambda(:,2)'*nor(:,3));
%   
%  
%   A(12,10) = QA(4,1)*Dlambda(:,3)'*nor(:,3);
%    
%   A(12,11) = QA(5,1)*Dlambda(:,3)'*nor(:,3);
%   A(12,13) = 2*QA(4,1)*Dlambda(:,2)'*nor(:,3) + 2*QA(5,1)*Dlambda(:,1)'*nor(:,3);
  %A(12,12) = -Dlambda(:,3)'*nor(:,3);
  %A(9,10)=QA(4,2)*(Dlambda(:,3)'*nor(:,3));
  
  %A(8,11)=QA(5,1)*(Dlambda(:,3)'*nor(:,3));
  %A(9,11)=QA(5,2)*(Dlambda(:,3)'*nor(:,3));
  
  %A(8,12)=0;
  %A(9,12)=0;
  
  %A(8,13)=2*QA(5,1)*(Dlambda(:,1)'*nor(:,3))+2*QA(4,1)*(Dlambda(:,2)'*nor(:,3));
  %A(9,13)=2*QA(5,2)*(Dlambda(:,1)'*nor(:,3))+2*QA(4,2)*(Dlambda(:,2)'*nor(:,3));
  
  %A(8,14)=3*QA(6,1)*(Dlambda(:,1)'*nor(:,3))+2*QA(7,1)*(Dlambda(:,2)'*nor(:,3));
  %A(9,14)=3*QA(6,2)*(Dlambda(:,1)'*nor(:,3))+2*QA(7,2)*(Dlambda(:,2)'*nor(:,3));
  
end

