function yyFilt = SGFilter(yy, k, n , s, varargin)
% Savitzky-Golay Filter
%
% Written By: Devon C. Hartlen, University of Waterloo, Jan 2022
%
% A smoothing methodology published by Savitzky and Golay (1964). This
% method fits a polynomial of order k to a window of n points around the 
% location of interest. For a 0-order polynomial, this filter is equivalent
% to a standard moving average filter. Additionally, the S-G filter can
% return the derivative over the same window of n points, up to a
% derivative order of k. 
%
% This implimentation of the S-G filter is does not crop or phase shift
% data. This is accomplished by special handling of beginning and end
% points. Performance drops off rapidly beyond fifth order polynomials. 
%
% INPUTS:
%   yy: 1D vector of points to be filtered. 
%   k: Order of the polynomial to be fit. Polynomial order must be less
%        than n-1, but should be much much less than n-1. A 0-order
%        polynomial is a standard moving average/boxcar filter. 
%   n: Number of points to average over (also called window size). Must be
%        an odd number
%   s: derivative order. If 0, the derivative is not computed. s cannot
%        exceed k. Requires dx to be defined to return unnormalized
%        values
%   dx: (optional). The physical spacing of points in the yy array.
%        Required for differentiation. If not provided, the normalized 
%        derivative is return and must be subsequently scaled by (dx^s)^-1  
%
% OUTPUT
%   yyFilt: Data smoothed or differentiated by the S-G filter

%% Input Handling
% Check if polynomial order exceeds window size
if k > (n-1)
    error(['Polynomial order k is equal to or greater than window size n', ...
        'k must not exceed n-1'])
end
% Check if window size is odd
if mod(n,2) ~= 1
    error('Window size n must be odd')
end
% Check that deriative order does not exceed polynomial order
if s > k
    error('Derivative order s must not exceed polynomial order k')
end

%% Perform Filtering
% Initalize array for filtered data
yyFilt = zeros(size(yy));
% Compute half-window size
m = (n-1)/2;
% Loop though all points in input array and perform filtering
for iPt = 1:length(yy)
    % This implimentation uses forward and backward fitting for the first
    % and last m points. This prevents cropping and phase shift. 
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

% If a deriviative was requested, check if a dx was provided. If so,
% pperform scaling. If not, warn the user and provide the normalized
% derivative
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

%% Ancillary Functions
% These functions are used to compute filter weights. Functions adapted
% from Gorry (1990), Anal. Chem., 60(5). 
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