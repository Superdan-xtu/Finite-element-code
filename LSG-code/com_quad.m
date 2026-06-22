function [ f ] = com_quad(alpha)
%COM_QUAD Summary of this function goes here
%   Detailed explanation goes here

n = size(alpha,1)-1;
for i = 1:size(alpha,1)
    if alpha(i)<0
        f=0;
        n=10000;
    end
end


switch n
case 1
f = factorial(alpha(1))*factorial(alpha(2))/factorial(sum(alpha)+n);
case 2
f = factorial(alpha(1))*factorial(alpha(2))*factorial(alpha(3))*factorial(n)/factorial(sum(alpha)+n);
case 3 
f = factorial(alpha(1))*factorial(alpha(2))*factorial(alpha(3))*factorial(alpha(4))*factorial(n)/factorial(sum(alpha)+n);
end
end

