function EQscale = eqScale(x,y)
% function to fit a exponential growth-saturation function to find the
% equilibrium scales of bedforms with... space or time... or whatever.

li = isnan(x) | isnan(y);

x(li) = [];
y(li) = []; 

x = x(:);
y = y(:);

p0 = [max(y)-abs(min(y)),1];

ft = fittype('a*(1-exp(-b*x))');
fo = fit(x,y,ft,'StartPoint', p0);

EQscale = fo.a;

