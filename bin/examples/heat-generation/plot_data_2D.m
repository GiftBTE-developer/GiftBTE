clc;clear all;close all
%%

load 001.dat;
Ta1=X001;
load 01.dat;
Ta2=X01;
load 05.dat
Ta3=X05;
Ta3(:,2)=Ta3(:,2)

semilogy(Ta1(:,1),Ta1(:,2),'b','LineWidth',2);
hold on


load 2e-7/TempLattice.dat;

C=1627100;
q=1E15;
k=150;
v=2677;
xx=0:0.01:1
S=q;
L=2e-7;
T11=TempLattice;
T11(:,4)=8*k*(T11(:,4))/S/L/L;
scatter(T11(:,2)/L,T11(:,4),80,'r');
hold on


plot(Ta2(:,1),Ta2(:,2),'b','LineWidth',2);
hold on


plot(Ta3(:,1),Ta3(:,2),'b','LineWidth',2);
hold on


load 1e-6/TempLattice.dat;

C=1627100;
q=1E15;
k=150;
v=2677;
xx=0:0.01:1
S=q;
L=1e-6;
T11=TempLattice;
T11(:,4)=8*k*(T11(:,4))/S/L/L;
scatter(T11(:,2)/L,T11(:,4),80,'r');
hold on


load 1e-5/TempLattice.dat;

C=1627100;
q=1E15;
k=150;
v=2677;
xx=0:0.01:1
S=q;
L=1e-5;
T11=TempLattice;
T11(:,4)=8*k*(T11(:,4))/S/L/L;
scatter(T11(:,2)/L,T11(:,4),80,'r');
hold on


box on

legend('Cao et al','Present','FontSize',24,'Location','best');
legend boxoff

xlabel('X*','FontSize',24);

ylabel('T*','FontSize',24);

set(gca,'FontSize',24);

text(0.5,2,'L=1000 nm','FontSize',24)

text(0.5,1,'L=100 nm','FontSize',24)

text(0.5,0.5,'L=10 nm','FontSize',24)
