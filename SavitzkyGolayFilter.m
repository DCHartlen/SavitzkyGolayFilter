fclose all;
close all;
clear;
clc;

cset = cbrewer2('set1',4)
colormap(cset)

k = 5; % polynomial order
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
yy = yyOg + noiseFact * (rand(nPts,1)-0.5*noiseFact);

figure();
subplot(1,2,1); hold on; 
plot(xx,yyOg, 'DisplayName','Ground Truth','color',cset(1,:))
plot(xx,yy, 'DisplayName','Added Noise','color',cset(2,:))
subplot(1,2,2); hold on;
plot(xx,yyOgDiff, 'DisplayName','Ground Truth','color',cset(1,:))
plot(xx(1:end-1)+0.5*dx,diff(yy)./dx, 'DisplayName','Added Noise','color',cset(2,:))

yyFilt = zeros(size(yy));
m = (winSize-1)/2;
for iPt = 1:length(yy)
    if iPt <= m
        t = iPt-m-1;
        weights = convWeights(k,s,t,winSize);
        yyFilt(iPt) = sum(yy(iPt-m-t:iPt+m-t).*weights);
    elseif iPt > length(yy)-m
        t = (iPt-length(yy))+m;
        weights = convWeights(k,s,t,winSize);
        yyFilt(iPt) = sum(yy(iPt-m-t:iPt+m-t).*weights);
    else
        weights = convWeights(k,s,0,winSize);
        yyFilt(iPt) = sum(yy(iPt-m:iPt+m).*weights);
    end
end
subplot(1,2,1); hold on; 
plot(xx,yyFilt, 'DisplayName', 'S-G Filter','color',cset(3,:))
legend()


yyFilt = zeros(size(yy));
m = (winSize-1)/2;
s=1;
for iPt = 1:length(yy)
    if iPt <= m
        t = iPt-m-1;
        weights = convWeights(k,s,t,winSize);
        yyFilt(iPt) = sum(yy(iPt-m-t:iPt+m-t).*weights);
    elseif iPt > length(yy)-m
        t = (iPt-length(yy))+m;
        weights = convWeights(k,s,t,winSize);
        yyFilt(iPt) = sum(yy(iPt-m-t:iPt+m-t).*weights);
    else
        weights = convWeights(k,s,0,winSize);
        yyFilt(iPt) = sum(yy(iPt-m:iPt+m).*weights);
    end
end
subplot(1,2,2); hold on; 
plot(xx,yyFilt./(dx^s), 'DisplayName','S-G Filt. Deriv','color',cset(3,:))
legend()

function weights = convWeights(k, s, t, winSize)
m = (winSize-1)/2;
weights = zeros(winSize,1);
for iIter = 1:winSize
    i = iIter-m-1;
    weights(iIter) = pointWeight(i, t, m, k, s);
end

end

function out = gramPolynomial(i, m, k, s)
if k > 0
    out = (4*k-2)/(k*(2*m-k+1))*(i*gramPolynomial(i, m, k-1, s) + s* ...
        gramPolynomial(i, m, k-1, s-1)) - ((k-1)*(2*m+k))/(k*(2*m-k+1))*...
        gramPolynomial(i, m, k-2, s);
elseif ((k==0) && (s==0))
    out = 1;
else
    out = 0;
end
end

function out = factors(a, b)
out = 1;
for iIter = (a-b+1):a
    out = out*iIter;
end
end

function weight = pointWeight(i, t, m, n, s)
weight = 0;
for k = 0:n
    weight = weight + (2*k+1)*(factors(2*m, k) / ...
        factors(2*m+k+1, k+1)) * gramPolynomial(i, m, k, 0) * ...
        gramPolynomial(t, m, k, s);
end
end