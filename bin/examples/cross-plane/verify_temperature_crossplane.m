clc;clear all;close all
%%
% compare temperature results from GiftBTE and analytical model for L=1e-8m
load 1e-8/TempLattice.dat; % load lattice temperature results calculated from GiftBTE
load analytic_201806/T1e-8.dat; % load results from analytical model
L=1e-8;
scatter(TempLattice(:,2)/L,TempLattice(:,4)-300,80,'r'); % non-dimensional x position versus Temperature rise
hold on

plot(T1e_8(:,1),T1e_8(:,2),'b','LineWidth',4); 
hold on

%%
% compare temperature results from GiftBTE and analytical model for L=1e-7m
load 1e-7/TempLattice.dat; % load lattice temperature results calculated from GiftBTE
load analytic_201806/T1e-7.dat; % load results from analytical model
L=1e-7;
scatter(TempLattice(:,2)/L,TempLattice(:,4)-300,80,'r'); % non-dimensional x position versus Temperature rise
hold on

plot(T1e_7(:,1),T1e_7(:,2),'b','LineWidth',4);
hold on

%%
% compare temperature results from GiftBTE and analytical model for L=1e-6m
load 1e-6/TempLattice.dat; % load lattice temperature results calculated from GiftBTE
load analytic_201806/T1e-6.dat; % load results from analytical model
L=1e-6;
scatter(TempLattice(:,2)/L,TempLattice(:,4)-300,80,'r'); % non-dimensional x position versus Temperature rise
hold on

plot(T1e_6(:,1),T1e_6(:,2),'b','LineWidth',4);
hold on

%%
box on

legend('GiftBTE','Analytical','FontSize',24);
legend boxoff

xlabel('X*','FontSize',24);

ylabel('T*','FontSize',24);

set(gca,'FontSize',24);

text(0.1,0.9,'L=1000 nm','FontSize',24)

text(0.3,0.7,'L=100 nm','FontSize',24)

text(0.1,0.5,'L=10 nm','FontSize',24)
