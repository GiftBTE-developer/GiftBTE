clc;clear all;close all
%%
load HeatFlux.dat; % load heat flux results from GiftBTE
figure
scatter3(HeatFlux(:,1),HeatFlux(:,2),HeatFlux(:,3),80,HeatFlux(:,4),'.'); % plot x-direction heat flux value
box on
colormap Turbo;
ch=colorbar;
set(get(ch,'title'),'string','unit:K');
set(gca,'FontSize',16);
xlabel('x(m)','FontSize',16);
ylabel('y(m)','FontSize',16);
zlabel('z(m)','FontSize',16);
title('x-HeatFlux (W/m^2)','FontSize',16);

%%
figure
scatter3(HeatFlux(:,1),HeatFlux(:,2),HeatFlux(:,3),80,HeatFlux(:,5),'.'); % plot y-direction heat flux value
box on
colormap Turbo;
ch=colorbar;
set(get(ch,'title'),'string','unit:K');
set(gca,'FontSize',16);
xlabel('x(m)','FontSize',16);
ylabel('y(m)','FontSize',16);
zlabel('z(m)','FontSize',16);
title('y-HeatFlux (W/m^2)','FontSize',16);

%%
figure
scatter3(HeatFlux(:,1),HeatFlux(:,2),HeatFlux(:,3),80,HeatFlux(:,6),'.'); % plot z-direction heat flux value
box on
colormap Turbo;
ch=colorbar;
set(get(ch,'title'),'string','unit:K');
set(gca,'FontSize',16);
xlabel('x(m)','FontSize',16);
ylabel('y(m)','FontSize',16);
zlabel('z(m)','FontSize',16);
title('z-HeatFlux (W/m^2)','FontSize',16);

%%
figure
% plot total heat flux value
scatter3(HeatFlux(:,1),HeatFlux(:,2),HeatFlux(:,3),80,(HeatFlux(:,6).^2+HeatFlux(:,5).^2+HeatFlux(:,4).^2).^0.5,'.'); 
box on
colormap Turbo;
ch=colorbar;
set(get(ch,'title'),'string','unit:K');
set(gca,'FontSize',16);
xlabel('x(m)','FontSize',16);
ylabel('y(m)','FontSize',16);
zlabel('z(m)','FontSize',16);
title('HeatFlux Magnitude (W/m^2)','FontSize',16);

