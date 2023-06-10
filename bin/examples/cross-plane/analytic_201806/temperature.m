temp=importdata("Tempcell.dat");
figure
plot(temp(:,1),temp(:,4));
hold on
x=0:0.01:1;
plot(x,T+300);