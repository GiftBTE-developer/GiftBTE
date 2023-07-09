clc;clear all;close all

%%
% load results from literature Ran, et al. (https://doi.org/10.1016/j.jcp.2019.108920)
load Tempcell_ref_20.dat; % L=20nm

load Tempcell_ref_50.dat; % L=50nm

load Tempcell_ref_220.dat; % L=220nm

%%
% compare temperature results from GiftBTE and literature Ran, et al.
plot(Tempcell_ref_20(:,1)/20e-9,Tempcell_ref_20(:,4),'b','LineWidth',2);
hold on
plot(Tempcell_ref_50(:,1)/50e-9,Tempcell_ref_50(:,4),'b','LineWidth',2);
hold on
h1=plot(Tempcell_ref_220(:,1)/220e-9,Tempcell_ref_220(:,4),'b','LineWidth',2);
hold on

load 20nm/TempLattice.dat; % GiftBTE L=20nm
scatter(TempLattice(:,2)/20e-9,TempLattice(:,4),80,'r'); % non-dimensional x position versus temperature
hold on

load 50nm/TempLattice.dat; % GiftBTE L=50nm
scatter(TempLattice(:,2)/50e-9,TempLattice(:,4),80,'r'); % non-dimensional x position versus temperature
hold on

load 220nm/TempLattice.dat; % GiftBTE L=220nm
h2=scatter(TempLattice(:,2)/220e-9,TempLattice(:,4),80,'r'); % non-dimensional x position versus temperature
hold on

%%
legend([h1 h2],'Ran et al','GiftBTE','FontSize',24,'Location','Southwest');
legend boxoff

xlabel('X*','FontSize',24);

ylabel('T*','FontSize',24);

set(gca,'FontSize',24);

text(0.5,0.9,'L=20 nm','FontSize',24)

text(0.5,0.8,'L=50 nm','FontSize',24)

text(0.5,0.6,'L=220 nm','FontSize',24)

box on
