close all; clear all; clc;
%% analytical result dirctly given
analytical5=load("analytical_result\kk5.txt"); % this file represent 2pikn=5
analytical1=load("analytical_result\kk1.txt"); % this file represent 2pikn=1
analytical2=load("analytical_result\kk2.txt"); % this file represent 2pikn=2
analytical05=load("analytical_result\kk05.txt"); % this file represent 2pikn=05
analytical025=load("analytical_result\kk025.txt"); % this file represent 2pikn=025
%% result from GiftBTE
num5=load("2pikn=5\TTG.dat"); % this file represent 2pikn=5
num1=load("2pikn=1\TTG.dat"); % this file represent 2pikn=1
num2=load("2pikn=2\TTG.dat"); % this file represent 2pikn=2
num05=load("2pikn=05\TTG.dat"); % this file represent 2pikn=05
num025=load("2pikn=025\TTG.dat"); % this file represent 2pikn=025

relaxationtime=3.99e-11; % relaxationtime from gray silicon band informantion (you can check it in "PHONON" file)

%% draw the temperature decay with time changing
plot(num5(:,1)/relaxationtime,num5(:,2),LineWidth=2,Color='r');hold on;%GitfBTE
plot(analytical5(:,1),analytical5(:,2),'--',LineWidth=2,Color='b');hold on;%analytical

plot(num2(:,1)/relaxationtime,num2(:,2),LineWidth=2,Color='r');hold on;%GitfBTE
plot(analytical2(:,1),analytical2(:,2),'--',LineWidth=2,Color='b');hold on;%analytical

plot(num1(:,1)/relaxationtime,num1(:,2),LineWidth=2,Color='r');hold on;%GitfBTE
plot(analytical1(:,1),analytical1(:,2),'--',LineWidth=2,Color='b');hold on;%analytical

plot(num05(:,1)/relaxationtime,num05(:,2),LineWidth=2,Color='r');hold on;%GitfBTE
plot(analytical05(:,1),analytical05(:,2),'--',LineWidth=2,Color='b');hold on;%analytical

plot(num025(:,1)/relaxationtime,num025(:,2),LineWidth=2,Color='r');hold on;%GitfBTE
plot(analytical025(:,1),analytical025(:,2),'--',LineWidth=2,Color='b');hold on;%analytical

xlim([0,10])
ylim([-0.2,1])

%%
box on

legend('GiftBTE','Analytical','FontSize',16);
legend boxoff

xlabel('t*','FontSize',24);

ylabel('T*','FontSize',24);

set(gca,'FontSize',24);

text(3.5,0.85,'2piKn=0.25','FontSize',16)

text(5.5,0.5,'2piKn=0.5','FontSize',16)

text(5,0.3,'2piKn=1','FontSize',16)

text(3,0.1,'2piKn=2','FontSize',16)
text(1.5,-0.1,'2piKn=5','FontSize',16)
