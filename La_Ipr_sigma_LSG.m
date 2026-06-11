function [ f1, f2, f3 ] = La_Ipr_sigma_LSG( X,La,mu )
%IPF Summary of this function goes here
%   Detailed explanation goes here
% f1 sigma_11, f2 sigma_22, f3 sigma_12
x = X(:,1);
y = X(:,2);
%f = sin(pi*x).*sin(pi*y);
% f1 =  3000.*pi.*cos(x.*pi).*sin(x.*pi).^2.*sin(y.*pi).^3;
% f2 =  1500.*pi.*cos(x.*pi).*sin(x.*pi).^2.*sin(y.*pi).^3 + 1500.*pi.*cos(y.*pi).*sin(x.*pi).^3.*sin(y.*pi).^2;
% f3 =  3000.*pi.*cos(y.*pi).*sin(x.*pi).^3.*sin(y.*pi).^2;
% f1 = f1/1000;
% f2 = f2/1000;
% f3 = f3/1000;
%f1 = 2.*mu.*((3.*(x - 1).^2.*(y - 3).^3)./La - 9.*pi.^2.*cos(x.*pi).*cos(y.*pi).*sin(x.*pi).^2.*sin(y.*pi).^2) + La.*((3.*(x - 1).^2.*(y - 3).^3)./La + (3.*(x - 1).^3.*(y - 3).^2)./La);
%f1 = La.*((y - 3).^2./La + ((2.*y - 6).*(x - 1))./La) + 2.*mu.*((y - 3).^2./La - 9.*pi.^2.*cos(x.*pi).*cos(y.*pi).*sin(x.*pi).^2.*sin(y.*pi).^2);
%f1 = 2.*mu.*((3.*(x - 1).^2.*(y - 1).^3)./(La + 1) - 9.*pi.^2.*cos(x.*pi).*cos(y.*pi).*sin(x.*pi).^2.*sin(y.*pi).^2) + La.*((3.*(x - 1).^2.*(y - 1).^3)./(La + 1) + (3.*(x - 1).^3.*(y - 1).^2)./(La + 1));
f1 =La.*((3.*x.^2.*y.^3.*(x - 1).^3.*(y - 1).^3)./(La + 1) + (3.*x.^3.*y.^2.*(x - 1).^3.*(y - 1).^3)./(La + 1) + (3.*x.^3.*y.^3.*(x - 1).^2.*(y - 1).^3)./(La + 1) + (3.*x.^3.*y.^3.*(x - 1).^3.*(y - 1).^2)./(La + 1)) + 2.*mu.*((3.*x.^2.*y.^3.*(x - 1).^3.*(y - 1).^3)./(La + 1) - 9.*pi.^2.*cos(x.*pi).*cos(y.*pi).*sin(x.*pi).^2.*sin(y.*pi).^2 + (3.*x.^3.*y.^3.*(x - 1).^2.*(y - 1).^3)./(La + 1));

f2 = La.*((3.*x.^2.*y.^3.*(x - 1).^3.*(y - 1).^3)./(La + 1) + (3.*x.^3.*y.^2.*(x - 1).^3.*(y - 1).^3)./(La + 1) + (3.*x.^3.*y.^3.*(x - 1).^2.*(y - 1).^3)./(La + 1) + (3.*x.^3.*y.^3.*(x - 1).^3.*(y - 1).^2)./(La + 1)) + 2.*mu.*(9.*pi.^2.*cos(x.*pi).*cos(y.*pi).*sin(x.*pi).^2.*sin(y.*pi).^2 + (3.*x.^3.*y.^2.*(x - 1).^3.*(y - 1).^3)./(La + 1) + (3.*x.^3.*y.^3.*(x - 1).^3.*(y - 1).^2)./(La + 1));

f3 =2.*mu.*(3.*pi.^2.*cos(x.*pi).^2.*sin(x.*pi).*sin(y.*pi).^3 - 3.*pi.^2.*cos(y.*pi).^2.*sin(x.*pi).^3.*sin(y.*pi) + (3.*x.^2.*y.^3.*(x - 1).^3.*(y - 1).^3)./(2.*(La + 1)) + (3.*x.^3.*y.^2.*(x - 1).^3.*(y - 1).^3)./(2.*(La + 1)) + (3.*x.^3.*y.^3.*(x - 1).^2.*(y - 1).^3)./(2.*(La + 1)) + (3.*x.^3.*y.^3.*(x - 1).^3.*(y - 1).^2)./(2.*(La + 1)));
%f1 = 2.*mu.*((pi.*cos(x.*pi).*sin(y.*pi))./(La + 1) - 9.*pi.^2.*cos(x.*pi).*cos(y.*pi).*sin(x.*pi).^2.*sin(y.*pi).^2) + La.*((pi.*cos(x.*pi).*sin(y.*pi))./(La + 1) + (pi.*cos(y.*pi).*sin(x.*pi))./(La + 1));
%f2 = 2.*mu.*((pi.*cos(y.*pi).*sin(x.*pi))./(La + 1) + 9.*pi.^2.*cos(x.*pi).*cos(y.*pi).*sin(x.*pi).^2.*sin(y.*pi).^2) + La.*((pi.*cos(x.*pi).*sin(y.*pi))./(La + 1) + (pi.*cos(y.*pi).*sin(x.*pi))./(La + 1));

%f3 = 2.*mu.*(3.*pi.^2.*cos(x.*pi).^2.*sin(x.*pi).*sin(y.*pi).^3 - 3.*pi.^2.*cos(y.*pi).^2.*sin(x.*pi).^3.*sin(y.*pi) + (pi.*cos(x.*pi).*sin(y.*pi))./(2.*(La + 1)) + (pi.*cos(y.*pi).*sin(x.*pi))./(2.*(La + 1)));

end

