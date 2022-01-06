function [yyFilt, A, B] = SaeJ211Filter(yy,cfc,fS)
% SAE J211-1 compliant filter
%
% Written By: Devon C. Hartlen, University of Waterloo, Jan 2022
%
% This is a four-pole, zero-phase, low-pass filter implimented per the SAE
% J211-1 standard.  
%
% INPUTS:
%   yy: 1D array of points to be filtered
%   cfc: cutoff frequency (Hz)
%   fS: sampling frequency of data (Hz)
% OUTPUTS:
%   filtered: Filtered 1D Array
%   A, B: filter constants (for reference only)

% Ensure Nyquist criteria is observed
if cfc > 0.5*fS
    error(['Nyquist criteria violated. Cutoff frequency must be less than' ...
        'half the sampling freqency'])
end

% Compute intermediate values
corner = 0.6*cfc;
dt = 1/fS;
wd = 2*pi*corner*2.0775;
wa = sin(wd*dt/2)/cos(wd*dt/2);

% Compute A coefficients
a0 = wa^2/(1+sqrt(2)*wa+wa^2);
a1 = 2*a0;
a2 = a0;

% Compute B coefficients
b0 = 1;
b1 = -2*(wa^2-1)/(1+sqrt(2)*wa+wa^2);
b2 = (-1+sqrt(2)*wa-wa^2)/(1+sqrt(2)*wa+wa^2);

A = [a0, a1, a2];
B = [b0, -b1, -b2];
% Function filtfilt performs efficient forward-backward convolution
% filtering.
yyFilt = filtfilt(A, B, yy);

end