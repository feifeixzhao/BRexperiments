function [medians, percentiles25, percentiles75] = calculateTemporalStats(runNumbers, individualDigitsArray)
    numSets = numel(runNumbers);
    medians = cell(1, numSets);
    percentiles25 = cell(1, numSets);
    percentiles75 = cell(1, numSets);

    for setIndex = 1:numSets
        runNumber = runNumbers(setIndex);
        individualDigits = individualDigitsArray{setIndex};

        digitsToProcess = unique([individualDigits]);

        numFiles = numel(digitsToProcess);
        medians{setIndex} = zeros(1, numFiles);
        percentiles25{setIndex} = zeros(1, numFiles);
        percentiles75{setIndex} = zeros(1, numFiles);

        folderName = sprintf('run%d', runNumber);

        for digitsIndex = 1:numFiles
            digits = digitsToProcess(digitsIndex);
            fileName = sprintf('%s_processedElevation_%03d.mat', folderName, digits);

            if exist(fullfile(folderName, fileName), 'file')
                fprintf('Processing file: %s\n', fileName);
                %load(fullfile(folderName, fileName));


                Zlid = scourspace(fullfile(folderName, fileName));
                Zlid(Zlid <= 0) = NaN;
          

                medians{setIndex}(digitsIndex) = nanmedian(Zlid(:));
                percentiles25{setIndex}(digitsIndex) = prctile(Zlid(:), 25, 'all'); % Calculate 25th percentile
                percentiles75{setIndex}(digitsIndex) = prctile(Zlid(:), 75, 'all'); % Calculate 75th percentile
            else
                fprintf('File not found: %s\n', fileName);
            end
        end
    end
end

