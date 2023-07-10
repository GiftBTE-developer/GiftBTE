clear all;clc;close all;
%%
load 1e-8/HeatFlux.dat; % load Heat Flux results calculated from GiftBTE for L=1e-8m
L=1e-8;
t=HeatFlux;
X1=0:0.01:1;
Y1=0:0.01:1;
[X1,Y1,Z]=griddata(t(:,1)/L,t(:,2)/L,t(:,4),X1',Y1,'linear');

figure
% contourf(X1,Y1,Z,0);
[C,h]=contourf(X1,Y1,Z,100); % plot x-direction Heat Flux contour
set(h,'LineColor','none');
colormap Turbo;
ch=colorbar;
set(get(ch,'title'),'string','W/m^2','Position',[40 -30]);
set(gca,'FontSize',20);
xlabel('X*','FontSize',20);
ylabel('Y*','FontSize',20);
title(['Heat Flux, L=',num2str(L),'m'],'FontSize',18);

%%
load 1e-7/HeatFlux.dat; % load Heat Flux results calculated from GiftBTE for L=1e-7m
L=1e-7;
t=HeatFlux;
X1=0:0.01:1;
Y1=0:0.01:1;
[X1,Y1,Z]=griddata(t(:,1)/L,t(:,2)/L,t(:,4),X1',Y1,'linear');

figure
% contourf(X1,Y1,Z,0);
[C,h]=contourf(X1,Y1,Z,100); % plot x-direction Heat Flux contour
set(h,'LineColor','none');
colormap Turbo;
ch=colorbar;
set(get(ch,'title'),'string','W/m^2','Position',[40 -30]);
set(gca,'FontSize',20);
xlabel('X*','FontSize',20);
ylabel('Y*','FontSize',20);
title(['Heat Flux, L=',num2str(L),'m'],'FontSize',18);

%%
load 1e-6/HeatFlux.dat; % load Heat Flux results calculated from GiftBTE for L=1e-6m
L=1e-6;
t=HeatFlux;
X1=0:0.01:1;
Y1=0:0.01:1;
[X1,Y1,Z]=griddata(t(:,1)/L,t(:,2)/L,t(:,4),X1',Y1,'linear');

figure
% contourf(X1,Y1,Z,0);
[C,h]=contourf(X1,Y1,Z,100); % plot x-direction Heat Flux contour
set(h,'LineColor','none');
colormap Turbo;
ch=colorbar;
set(get(ch,'title'),'string','W/m^2','Position',[40 -30]);
set(gca,'FontSize',20);
xlabel('X*','FontSize',20);
ylabel('Y*','FontSize',20);
title(['Heat Flux, L=',num2str(L),'m'],'FontSize',18);