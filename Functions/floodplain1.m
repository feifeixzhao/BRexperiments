function [Lid] = floodplain1(data_raw, x1, x2)
% calculates the average floodplain value of Ajay's braided river DEMs
% (isolating space)
% written by Feifei Zhao, 2021
%   Lid is the final average value, data_raw is the raw DEM data,
%   min_detrended is the reference point for y axis scaling
%   take data points 200 less than x1 and 200 more than x2
    data1=data_raw(x1-100:x1);
    data2=data_raw(x2:x2+100);
    data_raw=([data1;data2]);
    Lid = mean(data_raw);
    
end

