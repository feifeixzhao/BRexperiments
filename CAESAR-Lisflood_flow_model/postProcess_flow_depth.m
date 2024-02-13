%% plot colored data and histogram separately
run_number = 4; % Specify the run number
depth_time_number = 11; % Specify the modeled depth time number
elevation_time_number = 13; % Specify the elevation time number, +2 for Run2, +1 for others
M=postprocessmodel(run_number, depth_time_number, elevation_time_number);

percentiles = prctile(M.grid(:), [0, 25, 50, 75, 100]);

% Define a custom colormap with four colors on a red to blue gradient
cmap = [1 0 0;   % red  
        1 0.5 0.5; 
        0.5 0.5 1; 
        0 0 1];   %blue


% Apply the colormap to your data
colored_data = NaN(size(M.grid));
colored_data(M.grid < percentiles(2)) = 1;
colored_data(M.grid >= percentiles(4)) = 4;   
colored_data(M.grid >= percentiles(2) & M.grid < percentiles(3)) = 2; 
colored_data(M.grid >= percentiles(3) & M.grid < percentiles(4)) = 3; 
 

% Create the figure and plot the colored data
figure
imagesc(M.x, M.y, colored_data, 'alphadata',~isnan(M.grid))
xlim([15 30])
ylim([0.9 1.9])

colormap(cmap)
% Create a colorbar with custom labels
c = colorbar;
%c.TickLabels = cellstr(num2str(percentiles(1:5).'*1000)); % Set labels as the percentiles
c.Ticks = linspace(0, 5, 5); % Set the tick positions based on your data range

c.TickLabels = {'0', '25', '50', '75', '100'};

% Set the colorbar title
c.Label.String = 'percentile';

% plot corresponding histogram 
% Define the data and percentiles
data = M.grid(:).*1000; % Convert M.grid to a column vector
% Define the number of bins
num_bins = 30;

% Calculate the bin edges using the percentiles
percentiles = prctile(data, [0, 25, 50, 75, 100]);
bin_edges = linspace(percentiles(1), percentiles(end), num_bins + 1);

% Initialize the color array
cmap = [1 0 0; 1 0.5 0.5; 0.5 0.5 1; 0 0 1];
colors = zeros(num_bins, 3); % Initialize the colors array

% Calculate the 25th, 50th, and 75th percentiles of the entire data
data_percentiles = prctile(data, [25, 50, 75]);

% Assign colors to each bin based on global percentiles
for i = 1:num_bins
    if bin_edges(i+1) <= data_percentiles(1)
        colors(i, :) = cmap(1, :); % First color
    elseif bin_edges(i) >= data_percentiles(3)
        colors(i, :) = cmap(4, :); % Fourth color
    elseif bin_edges(i) >= data_percentiles(2) && bin_edges(i) < data_percentiles(3)
        colors(i, :) = cmap(3, :); % Third color
    elseif bin_edges(i) >= data_percentiles(1) && bin_edges(i) < data_percentiles(2)
        colors(i, :) = cmap(2, :); % Second color
    else
        colors(i, :) = cmap(1, :); % First color
    end
end

% Create the histogram using bar
figure;
[h, edges] = histcounts(data, num_bins); 
b = bar((edges(1:end-1)+edges(2:end))/2, h);
b.BarWidth=1; 
b.CData=colors;
b.FaceColor = 'flat';
xlim([0 15])

%% plot colored data and histogram together
run_numbers = [1, 2, 3, 4];

% Define the depth_time_number and elevation_time_number for each run
depth_time_numbers = [010, 005, 018, 009]; % Modify as needed
elevation_time_numbers = [11, 8, 20, 11]; % Modify as needed

% plotting 
figure;

% Create an empty cell array to store the handles to the colorbars
colorbar_handles = cell(1, length(run_numbers));

for i = 1:length(run_numbers)
    run_number = run_numbers(i);
    depth_time_number = depth_time_numbers(i);
    elevation_time_number = elevation_time_numbers(i);

    % Call the postprocessmodel function
    M = postprocessmodel(run_number, depth_time_number, elevation_time_number);

    % Plot the colored data
    subplot(length(run_numbers), 2, (i - 1) * 2 + 1);
    percentiles = prctile(M.grid(:), [0, 25, 50, 75, 100]);
    cmap = [1 0 0; 1 0.5 0.5; 0.5 0.5 1; 0 0 1];
    colored_data = NaN(size(M.grid));
    colored_data(M.grid < percentiles(2)) = 1;
    colored_data(M.grid >= percentiles(2) & M.grid < percentiles(3)) = 2;
    colored_data(M.grid >= percentiles(3) & M.grid < percentiles(4)) = 3;
    colored_data(M.grid >= percentiles(4)) = 4;
    imagesc(M.x, M.y, colored_data, 'alphadata', ~isnan(M.grid))
    xlim([15 30])
    ylim([0.7 2.1])
    colormap(cmap)
    
    % Add y-axis label to the first colored data plot
    if i == 3
        ylabel('Distance (m)');
    end

    if i == 4
        xlabel('Distance (m)');
    end


    % Plot the corresponding histogram
    subplot(length(run_numbers), 2, (i - 1) * 2 + 2);
    data = M.grid(:) * 1000;
    num_bins = 30;
    percentiles = prctile(data, [0, 25, 50, 75, 100]);
    bin_edges = linspace(percentiles(1), percentiles(end), num_bins + 1);
    cmap = [1 0 0; 1 0.5 0.5; 0.5 0.5 1; 0 0 1];
    colors = zeros(num_bins, 3);
    data_percentiles = prctile(data, [25, 50, 75]);
    for j = 1:num_bins
        if bin_edges(j + 1) <= data_percentiles(1)
            colors(j, :) = cmap(1, :);
        elseif bin_edges(j) >= data_percentiles(3)
            colors(j, :) = cmap(4, :);
        elseif bin_edges(j) >= data_percentiles(2) && bin_edges(j) < data_percentiles(3)
            colors(j, :) = cmap(3, :);
        elseif bin_edges(j) >= data_percentiles(1) && bin_edges(j) < data_percentiles(2)
            colors(j, :) = cmap(2, :);
        else
            colors(j, :) = cmap(1, :);
        end
    end
    [h, edges] = histcounts(data, num_bins);
    b = bar((edges(1:end - 1) + edges(2:end)) / 2, h);
    b.BarWidth = 1;
    b.CData = colors;
    b.FaceColor = 'flat';
    xlim([0 25])
    
    % Store the colorbar handle for the first colored data plot
    if i == 1
        colorbar_handles{i} = colorbar('Position', [0.13 0.95 0.275 0.015], 'Orientation', 'horizontal');
        colorbar_handles{i}.Ticks = 1:5;
        caxis([1,5])
        colorbar_handles{i}.TickLabels = {'0','25', '50', '75', '100'};
        colorbar_handles{i}.Label.String = 'percentile';
    end
end

% Add a single x-label for the last histogram
xlabel('h_{flow} (mm)', 'Interpreter', 'tex');

% Adjust subplot spacing
subplot(length(run_numbers), 2, 1);


