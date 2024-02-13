function [timeHours] = dimensionalize(dim_t, runNum)
    % Converts dimensionless time into experimental run hours
    % Predefined variables
    D_50 = 0.42 / 1000; % in meters
    g = 9.8; % gravity due to acceleration

    % Initialize variables
    Q = 0; % Flow rate in Liters
    S = 0; % Slope

    if runNum == 1
        Q = 0.25; % in Liters
        S = 0.01; % slope
    elseif runNum == 2
        Q = 0.12;
        S = 0.01;
    elseif runNum == 3
        Q = 0.25;
        S = 0.02;
    elseif runNum == 4
        Q = 0.12;
        S = 0.02; 
    else
        error('Invalid runNum. Only runNum 1 through 4 are supported.');
    end

    % Solve for timeHours using the rearranged equation
    timeHours = (dim_t * (Q * S^(2/3))) / (3600 * g * D_50^2 * 1000);
end
