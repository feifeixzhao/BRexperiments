function result = PBfitscour(Zlid_values)
    % Check if Zlid_values is provided
    if nargin == 0
        error('At least one input argument is required.');
    end
    
    % Analyze Zlid values
    if nargin >= 1
        % Remove negative, zero, and NaN values
        Zlid_values = Zlid_values(Zlid_values > 0 & ~isnan(Zlid_values));
        
        % Fit a gamma distribution
        pd = fitdist(Zlid_values, 'Gamma');
        scale_param = pd.ParameterValues(2);
        result = scale_param;
    end
    
end
