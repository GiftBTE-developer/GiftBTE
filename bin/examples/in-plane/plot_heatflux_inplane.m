clc;clear all;close all
%%
% compare x-direction heatflux results from GiftBTE and analytical model
load 1e-6/HeatFlux.dat; % GiftBTE L=1e-6m
load H1e-6.dat; % analytical model
L=1e-6;
scatter(HeatFlux(:,2)/L,HeatFlux(:,4),80,'r'); % non-dimensional y position versus x-direction heatflux
hold on
plot(H1e_6(:,1),H1e_6(:,2),'b','LineWidth',4);
hold on


load 1e-7/HeatFlux.dat; % GiftBTE L=1e-7m
load H1e-7.dat; % analytical model
L=1e-7;
scatter(HeatFlux(:,2)/L,HeatFlux(:,4),80,'r'); % non-dimensional y position versus x-direction heatflux
hold on
plot(H1e_7(:,1),H1e_7(:,2),'b','LineWidth',4);
hold on


load 1e-8/HeatFlux.dat; % GiftBTE L=1e-7m
load H1e-8.dat; % analytical model
L=1e-8;
scatter(HeatFlux(:,2)/L,HeatFlux(:,4),80,'r'); % non-dimensional y position versus x-direction heatflux
hold on
plot(H1e_8(:,1),H1e_8(:,2),'b','LineWidth',4);
hold on

%%
legend('GiftBTE','Analytical','FontSize',24,'Location','best');
legend boxoff

xlabel('Y*','FontSize',24);

ylabel('q_x (W/m^2)','FontSize',24);

set(gca,'FontSize',24);

text(0.3,10e8,'L=10 nm','FontSize',24)

text(0.3,6e8,'L=100 nm','FontSize',24)

text(0.3,2e8,'L=1000 nm','FontSize',24)
box on
