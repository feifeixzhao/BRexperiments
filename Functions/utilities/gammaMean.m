function M = gammaMean(x)
%mean of MLE gamma distribution fit
    x = x(~isnan(x));
    a = gamfit(x);
    [M,~] = gamstat(a(1),a(2));
    