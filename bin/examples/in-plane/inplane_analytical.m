clc;clear all;close all
load Input-dispersion-relation-fp-15.dat
data=Input_dispersion_relation_fp_15(:,:);
% load Input_CVT.dat
% load Input-dispersion-relation-fp.dat
% data=Input_dispersion_relation_fp(:,:);
% data=Input_CVT(:,:);
L=1e-8;
v=data(:,1);
tau=data(:,2);
c=data(:,3);
gradient=1;

Y1=0:0.01:1;
qbulk=sum(1/3*c.*v.*v.*tau);
for j=1:101
    Y=Y1(j);
    q(j)=0;
    for i=1:15
        fun=@(x) c(i)*tau(i).*v(i)^2.*(1-x.^2).*(2-exp(-Y./x/(v(i)*tau(i)/L))-exp(-(1-Y)./x/(v(i)*tau(i)/L)))
        q(j)=q(j)+1/4*integral(fun,0,1); 
    end
end

plot(Y1,q/L,'LineWidth',2);
hold on

result=[Y1',q'/L];
%load heatfluxfiel.dat
%scatter(heatfluxfiel(:,2),heatfluxfiel(:,3))
%hold on
% load 1e-9/heatfluxfiel_16.dat
% scatter(heatfluxfiel_16(:,2),heatfluxfiel_16(:,3))
% hold on
% load 1e-9/heatfluxfiel_32.dat
% scatter(heatfluxfiel_32(:,2),heatfluxfiel_32(:,3))
% hold on
% load 1e-9/heatfluxfiel_64.dat
% scatter(heatfluxfiel_64(:,2),heatfluxfiel_64(:,3))
% hold on
% load 1e-9/heatfluxfiel_128.dat
% scatter(heatfluxfiel_128(:,2),heatfluxfiel_128(:,3))
% hold on
% load 1e-9/heatfluxfiel_256.dat
% scatter(heatfluxfiel_256(:,2),heatfluxfiel_256(:,3))
% hold on

% load test-case/16*16/heatfluxfiel.dat
% scatter(heatfluxfiel(:,2),heatfluxfiel(:,3))
% hold on
% load test-case/32*32/heatfluxfiel.dat
% scatter(heatfluxfiel(:,2),heatfluxfiel(:,3))
% hold on
% load test-case/64*64/heatfluxfiel.dat
% scatter(heatfluxfiel(:,2),heatfluxfiel(:,3))
% hold on
% load test-case/128*128/heatfluxfiel.dat
% scatter(heatfluxfiel(:,2),heatfluxfiel(:,3))
% hold on
% load test-case/256*256/heatfluxfiel.dat
% scatter(heatfluxfiel(:,2),heatfluxfiel(:,3))
% hold on
% load inplane_chuang_6
% 
% scatter(inplane_chuang_6(:,1),inplane_chuang_6(:,2));

xlabel('{\ity^*}','fontsize',22);
ylabel('{\itq^*}','fontsize',22);

legend('Analytical','present 16*16','present 32*32','present 64*64','present 128*128','present 256*256');
legend boxoff

set(gca,'FontSize',20);