function [ phi1,phi2,phi3 ] = Cursefunction( node, elem, area, pxy )
%CURSEFUNCTION 此处显示有关此函数的摘要
%   此处显示详细说明
x = pxy(:,1);
y = pxy(:,2);
v1 = node(elem(:,1),:);
v2 = node(elem(:,2),:);
v3 = node(elem(:,3),:);
x1 = v1(:,1);
x2 = v2(:,1);
x3 = v3(:,1);
y1 = v1(:,2);
y2 = v2(:,2);
y3 = v3(:,2);
phi1 = ((x.*y2)./2 - (x2.*y)./2 - (x.*y3)./2 + (x3.*y)./2 + (x2.*y3)./2 - (x3.*y2)./2)./area;
phi2 = ((x1.*y)./2 - (x.*y1)./2 + (x.*y3)./2 - (x3.*y)./2 - (x1.*y3)./2 + (x3.*y1)./2)./area;
phi3 = ((x.*y1)./2 - (x1.*y)./2 - (x.*y2)./2 + (x2.*y)./2 + (x1.*y2)./2 - (x2.*y1)./2)./area;
phi1 = abs(phi1);
phi2 = abs(phi2);
phi3 = abs(phi3);
end

