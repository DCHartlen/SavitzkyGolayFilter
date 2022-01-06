function yyFilt = SGFilter(yy, k, n , s, varargin)


yyFilt = zeros(size(yy));
m = (n-1)/2;
for iPt = 1:length(yy)
    if iPt <= m
        t = iPt-m-1;
        weights = convWeights(k,s,t,n);
        yyFilt(iPt) = sum(yy(iPt-m-t:iPt+m-t).*weights);
    elseif iPt > length(yy)-m
        t = (iPt-length(yy))+m;
        weights = convWeights(k,s,t,n);
        yyFilt(iPt) = sum(yy(iPt-m-t:iPt+m-t).*weights);
    else
        weights = convWeights(k,s,0,n);
        yyFilt(iPt) = sum(yy(iPt-m:iPt+m).*weights);
    end
end

if s > 0
    if size(varargin,1) == 1
        dx = varargin{1};
        yyFilt = yyFilt./(dx^s);
    else
        warning(['Derivative requested but dx not provided.'...
            ' Normalized derivative returned'])
    end
end

end


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