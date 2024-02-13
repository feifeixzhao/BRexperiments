function medians = differenceDEM(runNumbers, individualDigitsArray)
    numSets = numel(runNumbers);
    medians = zeros(1, numSets);
    
    for setIndex = 1:numSets
        runNumber = runNumbers(setIndex);
        digits = individualDigitsArray{setIndex};
        
        % Load the first DEM
        folderName1 = sprintf('run%d', runNumber);
        fileName1 = sprintf('%s_processedElevation_%03d.mat', folderName1, digits(1));
        load(fullfile(folderName1, fileName1));
        dem1 = data.elevationMasked;
        
        % Load the second DEM
        folderName2 = sprintf('run%d', runNumber);
        fileName2 = sprintf('%s_processedElevation_%03d.mat', folderName2, digits(2));
        load(fullfile(folderName2, fileName2));
        dem2 = data.elevationMasked;
        
        % Calculate the difference between non-NaN values in dem1
        diff = dem2 - dem1;
        
        % Get non-NaN values in dem1
        nonNaNIndices = ~isnan(dem1);
        nonNaNDifferences = diff(nonNaNIndices);
        
        % Calculate the median of non-NaN differences
        medians(setIndex) = nanmean(nonNaNDifferences);
    end
end

