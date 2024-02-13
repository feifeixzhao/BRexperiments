function plotGammaAndEmpiricalPDFs(folderNames, fileNames, cmap, baseDirectory, useSubfolder)

    % Create a figure to hold all plots
    figure;

    % Number of folders
    numFolders = numel(folderNames);

    % Initialize legend entries for median plots
    legendEntries = cell(numFolders, 1);
    
    % Loop through the folders
    for i = 1:numFolders
        folderName = folderNames{i};
        fileName = fileNames{i};
        color = cmap(i, :);

        if useSubfolder
            % Construct the subfolder path
            subfolderPath = fullfile(folderName, 'braiding_phase', fileName);
        else
            % Use the main folder as the path
            subfolderPath = fullfile(folderName, fileName);
        end

        % Load Zclean data from the specified file
        filePath = fullfile(baseDirectory, subfolderPath);
        load(filePath);

        % Fit Zclean to a gamma distribution
        Pd = fitdist(Zclean, 'Gamma');

        % Create a range of values for the x-axis
        l = linspace(0, max(Zclean), 30);

        % Calculate the probability density function for the fitted gamma distribution
        gamma_pdf = pdf(Pd, l);

        %Zrnd = gamrnd(Pd.ParameterValues(1), Pd.ParameterValues(2), size(Zclean, 1), 1);
        %Zrnd(Zrnd < 0.001) = [];

        Z=Zclean;

        % Calculate the observed PDF (histogram)
        num_bins = 30;
        [heights, locations] = histcounts(Z, num_bins);
        width = (max(Z) - min(Z)) / num_bins;
        heights = heights / (length(Z) * width);

        % Create a bar plot for the histogram with the specified color and transparency
        h = bar(locations(1:end-1) * 1000, heights, 'hist');
        set(h, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', 0.3); % Color with 10% transparency

        hold on;

        % Plot the fitted gamma distribution in the specified color
        gammaPlot=plot(l * 1000, gamma_pdf, 'Color', color, 'LineWidth', 2);

        % Create a legend entry
        legendEntries{i} = gammaPlot;

    end

    % Add labels and a combined legend
    xlabel('Channel depth, h_{channel} (mm)');
    ylabel('Probability Density');
    disp(legendEntries)
    % Add legend
    legend([legendEntries{1}, legendEntries{2}, legendEntries{3}, legendEntries{4}], ...
        'Run 1, Q = 0.25 L/s, S = 0.01', 'Run 2, Q = 0.12 L/s, S = 0.01', 'Run 3, Q = 0.25 L/s, S = 0.02', 'Run 4, Q = 0.12 L/s, S = 0.02');


end


