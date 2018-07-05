function [meanphase,std,p] = kiri_circStats(phase)
% kiri_circStats is used for statistical analysis of phase lag data, uses
% circular statistics to calculate the circular mean, the circular standard
% deviation and finds an approximation of the p value for Rayleighs test
% for uniformity of distribution
%
% Kiri Pullar, masters thesis 2009

n=length(phase);
a=2*pi*phase;
c=sum(cos(a))/n;
s=sum(sin(a))/n;
ang=atan2(s,c);
meanang=mod(ang,2*pi);
meanphase=rad2deg(meanang)/360;

std = sqrt(-2*log((sqrt(c^2+s^2))))/360; 

R=n*sqrt(c^2+s^2);
z=R^2/n;

% approx p value of rayleighs R
p = exp(sqrt(1+(4*n)+4*(n^2-R^2))-(1+(2*n)));

end

