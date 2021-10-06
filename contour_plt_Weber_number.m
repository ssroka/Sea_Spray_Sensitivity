ccc
r = [50:10:2000]*1e-6;
u = 1:80;

s = 0.07;
rho = 1;
[U R] = meshgrid(u,r);
z = rho*U.^2.*R./s;
[c,h] = contourf(u,r,z,100);
set(h,'edgecolor','none')
hold on
contour(u,r,z,[1 1]*6,'r')
colorbar
g = gca;
g.ColorScale = 'log';
addpath /Users/ssroka/Documents/MATLAB/util/othercolor
%%
colormap(othercolor('Accent3'))
