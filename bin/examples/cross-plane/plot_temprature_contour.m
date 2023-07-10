clear all;clc;close all;
%%
load 1e-8/TempLattice.dat; % load temperature results calculated from GiftBTE for L=1e-8m
L=1e-8;
t=TempLattice;
X1=0:0.01:1;
Y1=0:0.01:1;
[X1,Y1,Z]=griddata(t(:,1)/L,t(:,2)/L,t(:,4),X1',Y1,'linear');

figure
% contourf(X1,Y1,Z,0);
[C,h]=contourf(X1,Y1,Z,100); % plot temperature contour
set(h,'LineColor','none');
colormap Turbo;
ch=colorbar;
set(get(ch,'title'),'string','unit:K');
set(gca,'FontSize',20);
xlabel('X*','FontSize',20);
ylabel('Y*','FontSize',20);
title(['Temperature, L=',num2str(L),'m'],'FontSize',18);

%%
load 1e-7/TempLattice.dat; % load temperature results calculated from GiftBTE for L=1e-7m
L=1e-7;
t=TempLattice;
X1=0:0.01:1;
Y1=0:0.01:1;
[X1,Y1,Z]=griddata(t(:,1)/L,t(:,2)/L,t(:,4),X1',Y1,'linear');

figure
% contourf(X1,Y1,Z,0);
[C,h]=contourf(X1,Y1,Z,100); % plot temperature contour
set(h,'LineColor','none');
colormap Turbo;
ch=colorbar;
set(get(ch,'title'),'string','unit:K');
set(gca,'FontSize',20);
title(['Temperature, ',num2str(L)]);
xlabel('X*','FontSize',20);
ylabel('Y*','FontSize',20);
title(['Temperature, L=',num2str(L),'m'],'FontSize',18);

%%
load 1e-6/TempLattice.dat; % load temperature results calculated from GiftBTE for L=1e-6m
L=1e-6;
t=TempLattice;
X1=0:0.01:1;
Y1=0:0.01:1;
[X1,Y1,Z]=griddata(t(:,1)/L,t(:,2)/L,t(:,4),X1',Y1,'linear');

figure
% contourf(X1,Y1,Z,0);
[C,h]=contourf(X1,Y1,Z,100); % plot temperature contour
set(h,'LineColor','none');
colormap Turbo;
ch=colorbar;
set(get(ch,'title'),'string','unit:K');
set(gca,'FontSize',20);
title(['Temperature, ',num2str(L)]);
xlabel('X*','FontSize',20);
ylabel('Y*','FontSize',20);
title(['Temperature, L=',num2str(L),'m'],'FontSize',18);
