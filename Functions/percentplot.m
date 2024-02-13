function percentplot(M)
    % Extract variables
    x = M.x;
    y = M.y;
    gridData = M.grid;

    % Calculate percentiles for coloring
    prctiles = prctile(gridData(:), [25, 50, 75]);
    
    % Create colormap based on percentiles
    cmap = [0 0 1; 0 1 0; 1 0 0]; % Blue, Green, Red

    % Normalize data to [0, 1] for colormap mapping
    normalizedData = (gridData - min(gridData(:))) / (max(gridData(:)) - min(gridData(:)));

    % Find which color each data point should be based on percentiles
    colorIndices = discretize(normalizedData, [0, prctiles, 1]);

    % Plot the DEM with colormap
    imagesc(x, y, gridData * 1000, 'alphadata', ~isnan(gridData));
    colormap(cmap);
    caxis([min(gridData(:)) * 1000, max(gridData(:)) * 1000]);
    colorbar;

    % Adjust colormap ticks and labels based on percentiles
    colorTicks = (prctiles - min(gridData(:))) * 1000;
    colorTickLabels = arrayfun(@(p) sprintf('%.2f', p), prctiles, 'UniformOutput', false);
    set(colorbar, 'Ticks', colorTicks, 'TickLabels', colorTickLabels);

    xlabel('X');
    ylabel('Y');
    title('DEM with Percentile-based Coloring');
end



