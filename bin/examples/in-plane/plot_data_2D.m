clc;clear all;close all
%%
load 2D-x/1e-6/HeatFlux.dat;
load H1e-6.dat;
L=1e-6;
scatter(HeatFlux(:,2)/L,HeatFlux(:,4),80,'r');
hold on


plot(H1e_6(:,1),H1e_6(:,2),'b','LineWidth',4);
hold on


load 2D-x/1e-7/HeatFlux.dat;
load H1e-7.dat;
L=1e-7;
scatter(HeatFlux(:,2)/L,HeatFlux(:,4),80,'r');
hold on



plot(H1e_7(:,1),H1e_7(:,2),'b','LineWidth',4);
hold on


load 2D-x/1e-8/HeatFlux.dat;
load H1e-8.dat;
L=1e-8;
scatter(HeatFlux(:,2)/L,HeatFlux(:,4),80,'r');
hold on



plot(H1e_8(:,1),H1e_8(:,2),'b','LineWidth',4);
hold on

legend('Present','Analytical','FontSize',24,'Location','best');
legend boxoff

xlabel('X*','FontSize',24);

ylabel('q (W/m^2)','FontSize',24);

set(gca,'FontSize',24);

text(0.3,10e8,'L=10 nm','FontSize',24)

text(0.3,6e8,'L=100 nm','FontSize',24)

text(0.3,2e8,'L=1000 nm','FontSize',24)
box on
