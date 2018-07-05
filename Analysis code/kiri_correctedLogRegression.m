function [loga, a,b, bCI, SE, rsquared,p] = kiri_correctedLogRegression(x,y)
% kiri_correctedLogRegression log transforms data, then examines correlation 
% and between variables. Values for a are corrected for log transformation 
% following Sprugel 1983. Values for log(a), a, b and its 95% confidence
% interval, the standard error associated with predictions, rsquared and
% p-values are returned.
%
% Kiri Pullar, masters thesis 2009
N=length(y);
logx=log(x);
logy=log(y);

ba = [logx,ones(size(logx))]\logy;
ypred = ba(1)*logx + ba(2);
SSE=sum((logy-ypred).^2);
SST=sum((logy-mean(logy)).^2);
SEE=sqrt((sum((logy-ypred).^2))/(N-2));
CF=exp((SEE)^2/2);
loga=ba(2);
a=exp(loga);
a=a*CF;
b=ba(1);

rsquared=1- SSE/SST;
tcrit=abs(tinv(0.025,N-2));
SE=sqrt((sum((logy-ypred).^2))/(N-2))/sqrt(sum((logx-mean(logx)).^2));
t=real(b/SE);
p=1-tcdf(t,N-2);
ME=tcrit*SE;
bCI=[b-ME b+ME]; %95% ci for slope