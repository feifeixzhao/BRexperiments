function cov = gammaCov(x)
%CoV for MLE gamma distribution fit
    x = x(~isnan(x));
    a = gamfit(x);
    [M,V] = gamstat(a(1),a(2));
    cov = sqrt(V)./M;