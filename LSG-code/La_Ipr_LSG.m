function [ f1, f2 ] = La_Ipr_LSG( X,La )
%IPF Summary of this function goes here
%   Detailed explanation goes here
x = X(:,1);
y = X(:,2);
%f = sin(pi*x).*sin(pi*y);
%f =  sin(pi*x).^2.*sin(pi*y).^2;
%f1 = 2.*mu.*((3.*(2.*x - 2).*(y - 3).^3)./La + 9.*pi.^3.*cos(y.*pi).*sin(x.*pi).^3.*sin(y.*pi).^2 - 18.*pi.^3.*cos(x.*pi).^2.*cos(y.*pi).*sin(x.*pi).*sin(y.*pi).^2) + 2.*mu.*((3.*(2.*y - 6).*(x - 1).^3)./(2.*La) + (9.*(x - 1).^2.*(y - 3).^2)./(2.*La) - 3.*pi.^3.*cos(y.*pi).^3.*sin(x.*pi).^3 + 6.*pi.^3.*cos(y.*pi).*sin(x.*pi).^3.*sin(y.*pi).^2 + 9.*pi.^3.*cos(x.*pi).^2.*cos(y.*pi).*sin(x.*pi).*sin(y.*pi).^2) + La.*((3.*(2.*x - 2).*(y - 3).^3)./La + (9.*(x - 1).^2.*(y - 3).^2)./La);
%f2 = 2.*mu.*((3.*(2.*y - 6).*(x - 1).^3)./La - 9.*pi.^3.*cos(x.*pi).*sin(x.*pi).^2.*sin(y.*pi).^3 + 18.*pi.^3.*cos(x.*pi).*cos(y.*pi).^2.*sin(x.*pi).^2.*sin(y.*pi)) + 2.*mu.*((3.*(2.*x - 2).*(y - 3).^3)./(2.*La) + (9.*(x - 1).^2.*(y - 3).^2)./(2.*La) + 3.*pi.^3.*cos(x.*pi).^3.*sin(y.*pi).^3 - 6.*pi.^3.*cos(x.*pi).*sin(x.*pi).^2.*sin(y.*pi).^3 - 9.*pi.^3.*cos(x.*pi).*cos(y.*pi).^2.*sin(x.*pi).^2.*sin(y.*pi)) + La.*((3.*(2.*y - 6).*(x - 1).^3)./La + (9.*(x - 1).^2.*(y - 3).^2)./La);
%f1 = ((x - 1).^3.*(y - 3).^3)./La - 3.*pi.*cos(y.*pi).*sin(x.*pi).^3.*sin(y.*pi).^2;
%f2 = ((x - 1).^3.*(y - 3).^3)./La + 3.*pi.*cos(x.*pi).*sin(x.*pi).^2.*sin(y.*pi).^3;
%f1 = ((x - 1).^3.*(y - 1).^3)./(La + 1) - 3.*pi.*cos(y.*pi).*sin(x.*pi).^3.*sin(y.*pi).^2;
%f2 = ((x - 1).^3.*(y - 1).^3)./(La + 1) + 3.*pi.*cos(x.*pi).*sin(x.*pi).^2.*sin(y.*pi).^3;
f1 = (x.^3.*y.^3.*(x - 1).^3.*(y - 1).^3)./(La + 1) - 3.*pi.*cos(y.*pi).*sin(x.*pi).^3.*sin(y.*pi).^2;
f2 = 3.*pi.*cos(x.*pi).*sin(x.*pi).^2.*sin(y.*pi).^3 + (x.^3.*y.^3.*(x - 1).^3.*(y - 1).^3)./(La + 1);
%f1 = (sin(x.*pi).*sin(y.*pi))./(La + 1) - 3.*pi.*cos(y.*pi).*sin(x.*pi).^3.*sin(y.*pi).^2;
%f2 = (sin(x.*pi).*sin(y.*pi))./(La + 1) + 3.*pi.*cos(x.*pi).*sin(x.*pi).^2.*sin(y.*pi).^3;
end

