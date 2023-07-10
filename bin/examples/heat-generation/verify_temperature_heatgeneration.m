clc;clear all;close all
%%
% load results from literature Cao, et al. (https://doi.org/10.1098/rspa.2015.0811)
load 001.dat; % load temperature results for L=1e-5m
Ta1=X001;
load 01.dat; % load temperature results for L=1e-6m
Ta2=X01;
load 05.dat  % load temperature results for L=2e-7m
Ta3=X05;

semilogy(Ta1(:,1),Ta1(:,2),'b','LineWidth',2);
hold on
plot(Ta2(:,1),Ta2(:,2),'b','LineWidth',2);
hold on
h1=plot(Ta3(:,1),Ta3(:,2),'b','LineWidth',2);
hold on

%%
% compare temperature results from GiftBTE and literature Cao, et al.
C=1627100;
q=1E15;
k=150;
v=2677;
xx=0:0.01:1
S=q;

load 2e-7/TempLattice.dat; % load temperature results calculated GiftBTE for L=2e-7m
L=2e-7;
T11=TempLattice;
T11(:,4)=8*k*(T11(:,4))/S/L/L;
scatter(T11(:,2)/L,T11(:,4),80,'r'); % non-dimensional x position versus temperature
hold on

load 1e-6/TempLattice.dat; % load temperature results calculated GiftBTE for L=1e-6m
L=1e-6;
T11=TempLattice;
T11(:,4)=8*k*(T11(:,4))/S/L/L;
scatter(T11(:,2)/L,T11(:,4),80,'r'); % non-dimensional x position versus temperature
hold on

load 1e-5/TempLattice.dat; % load temperature results calculated GiftBTE for L=1e-5m
L=1e-5;
T11=TempLattice;
T11(:,4)=8*k*(T11(:,4))/S/L/L;
h2=scatter(T11(:,2)/L,T11(:,4),80,'r'); % non-dimensional x position versus temperature
hold on

%%
box on

legend([h1 h2],'Cao et al','GiftBTE','Location','best');
% legend('Cao et al','GiftBTE','FontSize',24,'Location','best');
legend boxoff

xlabel('X*','FontSize',24);

ylabel('T*','FontSize',24);

set(gca,'FontSize',24);

text(0.5,2,'L=1000 nm','FontSize',24)

text(0.5,1,'L=100 nm','FontSize',24)

text(0.5,0.5,'L=10 nm','FontSize',24)
