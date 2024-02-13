function s = expStd(x)
%function to find the MLE fit of a exponential function
    x = x(~isnan(x));
    a = expfit(x);
    [~,V] = expstat(a);
    s = sqrt(V);
    