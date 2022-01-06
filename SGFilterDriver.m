fclose all;
close all;
clear;
clc;

cset = cbrewer2('set1',4);
colormap(cset)

k = 1; % polynomial order
winSize = 25; % window size (must be odd)
s = 0; % differentiation order

%% Fake experimental data
nPts = 500
xx = linspace(0,4,nPts)';
dx = (xx(end)-xx(1))./nPts;
yyOg = cos(pi.*xx);
yyOgDiff = -pi.*sin(pi.*xx);

rng(128493)
noiseFact = 0.5
yy = yyOg + noiseFact * (0.5-rand(nPts,1));

figure();
subplot(1,2,1); hold on; 
plot(xx,yyOg, 'DisplayName','Ground Truth','color',cset(1,:))
plot(xx,yy, 'DisplayName','Added Noise','color',cset(2,:))
subplot(1,2,2); hold on;
plot(xx,yyOgDiff, 'DisplayName','Ground Truth','color',cset(1,:))
plot(xx(1:end-1)+0.5*dx,diff(yy)./dx, 'DisplayName','Added Noise','color',cset(2,:))

yyFilt = SGFilter(yy,k, winSize,s);
subplot(1,2,1)
plot(xx,yyFilt, 'DisplayName','S-G Filter','color',cset(3,:))
legend()

s = 1;
yyDerivFilt = SGFilter(yy, k, winSize, s, dx);
subplot(1,2,2);
plot(xx,yyDerivFilt, 'DisplayName','S-G Filter Deriv.','color',cset(3,:))
