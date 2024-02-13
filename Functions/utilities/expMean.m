function s = expMean(x)
%function to find the MLE fit of an exponential distribution
    x = x(~isnan(x));
    a = expfit(x);
    [M,~] = expstat(a);
    s = M;
    