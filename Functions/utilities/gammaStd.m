function s = gammaStd(x)
% standard deviation from MLE fit of gamma distribution
    x = x(~isnan(x));
    a = gamfit(x);
    [~,V] = gamstat(a(1),a(2));
    s = sqrt(V);
    