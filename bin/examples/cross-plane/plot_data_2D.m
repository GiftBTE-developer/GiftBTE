clc;clear all;close all
%%
load 2D-y/1e-8/Tempcell.dat;
load analytic_201806/T1e-8.dat;
L=1e-8;
scatter(Tempcell(:,2)/L,Tempcell(:,4),80,'r');
hold on

plot(T1e_8(:,1),T1e_8(:,2),'b','LineWidth',4);
hold on

load 2D-y/1e-7/Tempcell.dat;
load analytic_201806/T1e-7.dat;
L=1e-7;
scatter(Tempcell(:,2)/L,Tempcell(:,4),80,'r');
hold on

plot(T1e_7(:,1),T1e_7(:,2),'b','LineWidth',4);
hold on


load 2D-y/1e-6/Tempcell.dat;
load analytic_201806/T1e-6.dat;
L=1e-6;
scatter(Tempcell(:,2)/L,Tempcell(:,4),80,'r');
hold on

plot(T1e_6(:,1),T1e_6(:,2),'b','LineWidth',4);
hold on

box on

legend('Analytical','Present','FontSize',24);
legend boxoff

xlabel('X*','FontSize',24);

ylabel('T*','FontSize',24);

set(gca,'FontSize',24);

text(0.1,0.9,'L=1000 nm','FontSize',24)

text(0.3,0.7,'L=100 nm','FontSize',24)

text(0.1,0.5,'L=10 nm','FontSize',24)
