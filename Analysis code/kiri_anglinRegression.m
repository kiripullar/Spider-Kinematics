function [rsquared,p] = kiri_anglinRegression(phase,x)
% kiri_anglinRegression examines angular-linear correlation and returns the
% rsquared value and p-value associated with this correlation.
%
% Kiri Pullar, masters thesis 2009

a=2*pi*phase;
n = length(a);

% correlation between x and cosine of the angle, correlation between  x  and
% the sine of the angle and correlation between the sine and cosine of the
% angle
rxs = corr(x,sin(a));
rxc = corr(x,cos(a));
rcs = corr(sin(a),cos(a));

% compute angular-linear correlation
rsquared = sqrt((rxc^2 + rxs^2 - 2*rxc*rxs*rcs)/(1-rcs^2));

% compute pvalue
p = chi2pdf(n*rsquared^2,2);