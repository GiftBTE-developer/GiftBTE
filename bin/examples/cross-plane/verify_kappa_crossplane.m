clc;clear all;close all
%%
% calculate effective thermal conductivity using results from GiftBTE and compare with analytical model
load analytic_201806/TC.dat; % load effective thermal conductivity results from analytical model
load 1e-8/HeatFlux.dat; % load heat flux results calculated from GiftBTE for L=1e-8m

L(1)=1e-8;
k(1)=sum(HeatFlux(:,5))/10000*L(1); % calculated effective thermal conductivity using heat flux calculated from GiftBTE

load 1e-7/HeatFlux.dat; % load heat flux results calculated from GiftBTE for L=1e-7m

L(2)=1e-7;
k(2)=sum(HeatFlux(:,5))/10000*L(2); % calculated effective thermal conductivity using heat flux calculated from GiftBTE

load 1e-6/HeatFlux.dat; % load heat flux results calculated from GiftBTE for L=1e-6m

L(3)=1e-6;
k(3)=sum(HeatFlux(:,5))/10000*L(3); % calculated effective thermal conductivity using heat flux calculated from GiftBTE

semilogx(TC(:,1),TC(:,2),'b','LineWidth',2); 
hold on;
scatter(L,k,200,'r^') % length L versus GiftBTE calculated effective thermal conductivity

%%
legend('Analytical','GiftBTE','FontSize',24,'Location','northwest');
legend boxoff

xlabel('L (m)','FontSize',24);

ylabel('\kappa (W/mK)','FontSize',24);

set(gca,'FontSize',24);

