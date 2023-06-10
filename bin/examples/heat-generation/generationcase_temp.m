clear all; close all;clc;


load 001.dat;
Ta1=X001;
load 01.dat;
Ta2=X01;
load 05.dat
Ta3=X05;
Ta3(:,2)=Ta3(:,2)

C=1627100;
q=1E15;
k=150;
v=2677;
xx=0:0.01:1
S=q;
L=2e-7;

load Tempcell.dat;
T11=Tempcell;
T11(:,4)=8*k*(T11(:,4))/S/L/L;

semilogy(Ta1(:,1),Ta1(:,2),'r','LineWidth',2);
hold on
scatter(T11(:,1)/L,T11(:,4),'MarkerEdgecolor','b');
hold on

plot(Ta2(:,1),Ta2(:,2),'r','LineWidth',2);
hold on


plot(Ta3(:,1),Ta3(:,2),'r','LineWidth',2);
hold on


% load testforsquare/internalheat_gray1/Tempcell.dat;
% T11=Tempcell;
% T11(:,3)=8*k*(T11(:,3))/S/1e-4/1e-4;
% load testforsquare_1200/internalheat_gray1/Tempcell.dat;
% T21=Tempcell;
% T21(:,3)=8*k*(T21(:,3))/S/1e-4/1e-4;
% load testforsquare_225/internalheat_gray1/Tempcell.dat;
% T31=Tempcell;
% T31(:,3)=8*k*(T31(:,3))/S/1e-4/1e-4;
% scatter(T11(:,1),T11(:,3),'MarkerEdgecolor','b');
% hold on
% scatter(T21(:,1),T21(:,3),'MarkerEdgecolor','k');
% hold on
% scatter(T31(:,1),T31(:,3),'MarkerEdgecolor','g');
% hold on


legend('Cao, 2016','Present');
legend boxoff

set(gca,'FontSize',20);