function [Zlid, xVec, yVec] = scourspace(matFilePath)
    % Load the data from the specified .mat file
    load(matFilePath);

    % Initialize variables to store the results
    numColumns = size(data.elevationMasked, 2);
    columnAverages = NaN(1, numColumns);

    % Initialize a new matrix to store masked raw data
    maskedRawData = NaN(size(data.elevationRaw));

    % Iterate through each column
    for col = 1:numColumns
        maskedColumn = data.elevationMasked(:, col);
        rawColumn = data.elevationRaw(:, col);
        
        % Find the indices of the first and last non-NaN values in the masked column
        firstNonNaNIndex = find(~isnan(maskedColumn), 1, 'first');
        lastNonNaNIndex = find(~isnan(maskedColumn), 1, 'last');
        
        % Check if the column has any valid data
        if ~isempty(firstNonNaNIndex) && ~isempty(lastNonNaNIndex)
            % Compute the average for the specified range in the raw column
            columnAverages(col) = nanmean([rawColumn(1:firstNonNaNIndex-1); rawColumn(lastNonNaNIndex+1:end)]);
            
            % Create a new column in maskedRawData that contains only the values between firstNonNaNIndex and lastNonNaNIndex
            maskedRawData(:, col) = rawColumn;
            maskedRawData(1:firstNonNaNIndex-1, col) = NaN;
            maskedRawData(lastNonNaNIndex+1:end, col) = NaN;
        end
    end

    % Find the indices of non-NaN values in columnAverages
    validIndices = ~isnan(columnAverages);
    numValidColumns = sum(validIndices);

    % Interpolate the columnAverages for valid columns only
    interpolatedValues = NaN(1, numColumns);
    if numValidColumns > 1
        firstValidIndex = find(validIndices, 1, 'first');
        lastValidIndex = find(validIndices, 1, 'last');
        interpolatedValues(firstValidIndex:lastValidIndex) = linspace(columnAverages(firstValidIndex), columnAverages(lastValidIndex), numValidColumns);
    end

    % Subtract the interpolated values from each row of data.elevationMasked
    interpolatedMatrix = repmat(interpolatedValues, size(data.elevationMasked, 1), 1);
    Zlid = interpolatedMatrix - maskedRawData;

    % Output xVec and yVec
    xVec = data.xVec;
    yVec = data.yVec;
end


