clc;clear all;close all

%%

load Tempcell_ref_20.dat;

load Tempcell_ref_50.dat;

load Tempcell_ref_220.dat;

plot(Tempcell_ref_20(:,1)/20e-9,Tempcell_ref_20(:,4),'b','LineWidth',2);
hold on

load 20nm/TempLattice.dat;

scatter(Tempcell(:,2)/20e-9,Tempcell(:,4),80,'r');
hold on

plot(Tempcell_ref_50(:,1)/50e-9,Tempcell_ref_50(:,4),'b','LineWidth',2);
hold on

plot(Tempcell_ref_220(:,1)/220e-9,Tempcell_ref_220(:,4),'b','LineWidth',2);
hold on



load 50nm/TempLattice.dat;

scatter(Tempcell(:,2)/50e-9,Tempcell(:,4),80,'r');
hold on

load 220nm/TempLattice.dat;

scatter(Tempcell(:,2)/220e-9,Tempcell(:,4),80,'r');
hold on

legend('Ran et al','Present','FontSize',24,'Location','best');
legend boxoff

xlabel('X*','FontSize',24);

ylabel('T*','FontSize',24);

set(gca,'FontSize',24);

text(0.5,0.9,'L=20 nm','FontSize',24)

text(0.5,0.8,'L=50 nm','FontSize',24)

text(0.5,0.6,'L=220 nm','FontSize',24)

box on
