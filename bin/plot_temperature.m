clc;clear all;close all
%%
load TempLattice.dat; % load temperature results from GiftBTE
figure
scatter3(TempLattice(:,1),TempLattice(:,2),TempLattice(:,3),80,TempLattice(:,4),'.'); % plot temperature
box on
colormap Turbo;
ch=colorbar;
set(get(ch,'title'),'string','unit:K');
set(gca,'FontSize',16);
xlabel('x(m)','FontSize',16);
ylabel('y(m)','FontSize',16);
zlabel('z(m)','FontSize',16);
title('Temperature','FontSize',16);

