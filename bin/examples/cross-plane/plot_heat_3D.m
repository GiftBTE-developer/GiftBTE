clc;clear all;close all

load analytic_201806/TC.dat;
load 3D-z/1e-8/HeatFlux.dat;

L(1)=1e-8;
k(1)=sum(HeatFlux(:,6))/8000*L(1);

load 3D-z/1e-7/HeatFlux.dat;

L(2)=1e-7;
k(2)=sum(HeatFlux(:,6))/8000*L(2);

load 3D-z/1e-6/HeatFlux.dat;

L(3)=1e-6;
k(3)=sum(HeatFlux(:,6))/8000*L(3);

semilogx(TC(:,1),TC(:,2),'b','LineWidth',2);
hold on;
scatter(L,k,200,'r^')

legend('Analytical','Present','FontSize',24,'Location','northwest');
legend boxoff

xlabel('L (m)','FontSize',24);

ylabel('\kappa (W/m-K)','FontSize',24);

set(gca,'FontSize',24);

