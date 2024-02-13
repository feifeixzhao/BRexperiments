% script is to code figures in prep for JSR paper (Zhao et al., 2024). 

%% Figure 1: plot DEMs from all four runs
run_numbers = [1, 2, 1, 3, 1, 4];

% Define the last three digits of the file names
file_numbers = [002,012, 007, 021, 032, 020];

fig=figure(1)
for i = 1:length(run_numbers)
    run_number = run_numbers(i);
    
    % Load the corresponding Zlid data based on the last three digits
    % of the file name
    file_number = file_numbers(i);
    filename = sprintf('/Users/Feifei/Box/BR_experiments/braidingExperimentsTopo_v2/run%d/run%d_processedElevation_%03d.mat', run_number, run_number, file_number);
    
    % Load the data as a structure
    load(filename); 
    
    % Create subplots for each run
    subplot(3, 2, i); 
    % Modify the y values
    yVec_median = median(data.yVec);
    yVec_new = data.yVec - yVec_median;
    imagesc(data.xVec(9500:12050), yVec_new(252:1100), (data.elevationMaskedDetrended(252:1100, 9500:12050)) * 1000, 'alphadata', ~isnan(data.elevationMaskedDetrended(252:1100, 9500:12050)));
    colormap(jet)
    caxis([-5 5])
end
h = axes(fig,'visible','off'); 
caxis([-5 5])
c = colorbar('north')
set(h,'Position',[0.35 0.95 0.3 0.062]);
%c = colorbar(h,'Location', 'Northoutside');
%c.Label.String='Detrended Elevation (mm)'


%% Figure 3: fit scour depths to gamma
baseDirectory = '/Users/Feifei/Box/BR_experiments/braidingExperimentsTopo_v2';

folderNames = {'run1_data', 'run2_data', 'run3_data', 'run4_data'};
fileNames = {'Zclean.mat', 'Zclean.mat', 'Zclean.mat', 'Zclean.mat'};
cmap = [
    0 0 0;        % Black
    0.4 0.4 0.4;  % Light Gray
    1 0 0;        % Red
    1 0.6 0.6;    % Light Red
];
useSubfolder = false;  % Set to true to access Zclean in the "braiding_phase" subfolder

plotGammaAndEmpiricalPDFs(folderNames, fileNames, cmap, baseDirectory, useSubfolder);
xlim([0 30])

%% goodness of fit test 
Z=Z_braiding.*1000; % convert scours to mm 
bins=30; 

l=linspace(0, 48, bins);
pd=fitdist(Z, 'gamma'); 
Z_pred=pdf(pd, l); 

% Normalize Z_pred so that the sum is equal to 1
Z_pred = Z_pred / sum(Z_pred);
binCenters = [l(1:end) - 0.5 * l(2), l(end)]; % add another value at the very end

% Calculate the observed histogram
Z_obs = histcounts(Z, binCenters, 'Normalization', 'probability');

% calculate chi squared stat
chi_2=sum((Z_obs-Z_pred).^2/Z_pred);

% Degrees of freedom
df = numel(l) - 1;

% P-value from chi-squared distribution
pValue = 1 - chi2cdf(chi_2, df);

% Display results
fprintf('Chi-squared test results:\n');
fprintf('Chi2 statistic = %.4f\n', chi_2);
fprintf('Degrees of freedom = %d\n', df);
fprintf('P-value = %.4f\n', pValue);

% Determine whether to reject or accept the null hypothesis
alpha = 0.05; % Set your significance level
if pValue < alpha
    fprintf('Reject the null hypothesis. The distributions are different.\n');
else
    fprintf('Accept the null hypothesis. The distributions are similar.\n');
end


%% Figure 4: plot modeled flow depth distributions 
% Define base folder and file name

% Define base folder
baseFolder = '/Users/feife/Box/BR_experiments/caesar_lisflood_2022/new_flow_depths';

% Search for folders that start with "run"
folderPattern = 'run*';
folderList = dir(fullfile(baseFolder, folderPattern));

% Define the number of bins for the histograms
numBins = 30;

% Initialize a figure
figure;

% Define the colors for the histograms
cmap = [
    0 0 0;        % Black
    0.5 0.5 0.5;  % Light Gray
    1 0 0;        % Red
    1 0.6 0.6;    % Light Red
];

% Loop through each folder
for folderIdx = 1:numel(folderList)
    folderName = folderList(folderIdx).name;
    filePath = fullfile(baseFolder, folderName, 'new_flow_depths.mat');
    
    % Load data from the .mat file
    data = load(filePath);
    
    % Combine M1 to M5 into a single variable
    combinedData = [data.M1(:); data.M2(:); data.M3(:); data.M4(:); data.M5(:)].*1000;
    
    % Define the number of bins for the histogram
    numBins = 30;

    % Plot relative frequency histogram for the combined data
    histogram(combinedData, numBins, 'Normalization', 'probability', 'FaceColor', cmap(folderIdx, :), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    
    % Adjust the hold state so that multiple histograms are overlaid
    hold on;
end

% Add labels and a legend
xlim([0 30])
xlabel('Modeled flow depths, h_{flow} (mm)');
ylabel('Frequency');
legend({'Run 1, Q = 0.25 L/s, S = 0.01', 'Run 2, Q = 0.12 L/s, S = 0.01', 'Run 3, Q = 0.25 L/s, S = 0.02', 'Run 4, Q = 0.12 L/s, S = 0.02'});

% Adjust the figure layout
set(gcf, 'Position', [100, 100, 800, 600]);

%% Figure 6: plot scour depth distribution 
% load and calculate scour depth values
% Initialize the figure
figure;
run_numbers = [1, 2, 3, 4];

% Create an empty cell array to store the handles to the colorbars
colorbar_handles = cell(1, length(run_numbers));

% Define the last three digits of the file names
file_numbers = [011, 008, 020, 010]; % Modify as needed
for i = 1:length(run_numbers)
    run_number = run_numbers(i);
    
     % Load the corresponding Zlid data based on the last three digits
    % of the file name
    file_number = file_numbers(i);
    filename = sprintf('/Users/feife/Box/BR_experiments/braidingExperimentsTopo_v2/run%d/run%d_processedElevation_%03d.mat', run_number, run_number, file_number);
    
    % Process Zlid data using the scourspace function
    [Zlid, xVec, yVec] = scourspace(filename); % Pass the file path as an argument
    Zlid(Zlid <= 0) = NaN;

    % Plot the colored data
    subplot(length(run_numbers), 2, (i - 1) * 2 + 1);
    percentiles = prctile(Zlid(:), [0, 25, 50, 75, 100]);
    cmap = [1 0 0; 1 0.5 0.5; 0.5 0.5 1; 0 0 1];
    colored_data = NaN(size(Zlid));
    colored_data(Zlid < percentiles(2)) = 1;
    colored_data(Zlid >= percentiles(2) & Zlid < percentiles(3)) = 2;
    colored_data(Zlid >= percentiles(3) & Zlid < percentiles(4)) = 3;
    colored_data(Zlid >= percentiles(4)) = 4;
    imagesc(xVec, yVec, colored_data, 'alphadata', ~isnan(Zlid))
    xlim([15.8 30.8])
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
    data = Zlid(:) * 1000;
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
        caxis([1, 5])
        colorbar_handles{i}.TickLabels = {'0', '25', '50', '75', '100'};
        colorbar_handles{i}.Label.String = 'percentile';
    end
end

% Add a single x-label for the last histogram
xlabel('h_{channel} (mm)', 'Interpreter', 'tex');

% Adjust subplot spacing
subplot(length(run_numbers), 2, 1);

%% Figure 5: plot scour depths and flow depths through time 
% load flow depths 
baseFolder = '/Users/feife/Box/BR_experiments/caesar_lisflood_2022/new_flow_depths';
fileName = 'new_flow_depths.mat'; % Assuming the file name is the same for all folders

% Search for folders that start with "run"
folderPattern = 'run*';
folderList = dir(fullfile(baseFolder, folderPattern));
numFolders = numel(folderList);

% Initialize arrays to store results
allMedians = [];
allPercentiles25 = [];
allPercentiles75 = [];

% Loop through each folder
for folderIdx = 1:numFolders
    folderName = folderList(folderIdx).name;
    filePath = fullfile(baseFolder, folderName, fileName);
    
    % Load data from the .mat file
    data = load(filePath);
    
    % Calculate median, 25th and 75th percentiles for each variable
    variables = fieldnames(data);
    medians = [];
    percentiles25 = [];
    percentiles75 = [];
    
    for varIdx = 1:numel(variables)
        variableName = variables{varIdx};
        variableData = data.(variableName);
        
        medianValue = median(variableData, 'all');
        percentile25Value = prctile(variableData, 25, 'all'); 
        percentile75Value = prctile(variableData, 75, 'all');
        
        medians(varIdx) = medianValue;
        percentiles25(varIdx) = percentile25Value;
        percentiles75(varIdx) = percentile75Value;
    end
    
    % Store results for this folder
    allMedians(folderIdx, :) = medians;
    allPercentiles25(folderIdx, :) = percentiles25;
    allPercentiles75(folderIdx, :) = percentiles75;
end

% Define dimensionless time
dimt = [0, 2000, 4000, 6000, 8000, 10000];

% Define custom colors
customColors = [
    0 0 0;        % Black
    0.5 0.5 0.5;  % Light Gray
    1 0 0;        % Red
    1 0.6 0.6;    % Light Red
];

% Create a new figure
figure;

% Initialize legend entries for median plots
legendEntries = cell(numFolders, 1);

% Loop through each row in allMedians
for rowIdx = 1:size(allMedians, 1)
    medianData = allMedians(rowIdx, :);
    percentile25Data = allPercentiles25(rowIdx, :);
    percentile75Data = allPercentiles75(rowIdx, :);
    color = customColors(rowIdx, :);
    
    % Plot the area between 25th and 75th percentiles with matching dashed lines
    fill([dimt, fliplr(dimt)], [percentile25Data.*1000, fliplr(percentile75Data.*1000)], color, 'FaceAlpha', 0.15, 'EdgeColor', color, 'LineStyle', '--');
    hold on
    
    % Plot median data with dashed lines and circles (adjust LineWidth here)
    medianPlot = plot(dimt, medianData.*1000, 'Color', color, 'LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor', color, 'LineWidth', 2);
    
    % Create a legend entry for the median plot
    legendEntries{rowIdx} = medianPlot;
end

% Add legend
legend([legendEntries{1}, legendEntries{2}, legendEntries{3}, legendEntries{4}], ...
  'Run 1, Q = 0.25 L/s, S = 0.01', 'Run 2, Q = 0.12 L/s, S = 0.01', 'Run 3, Q = 0.25 L/s, S = 0.02', 'Run 4, Q = 0.12 L/s, S = 0.02');


% Add labels and title
xlabel('Dimensionless time, {\it t*}');
ylabel('Channel depth, {\it h_{flow}} (mm)')

% Set x-axis tick values
xticks(xticks);
xticklabels(arrayfun(@(x) sprintf('%g', x), xticks/1000, 'UniformOutput', false));

% Remove top border
ax = gca; % Get the current axes
ax.Box = 'on'; 

%% Figure 5 plot scour depths over time 
runNumbers = [1, 2, 3, 4];  % Array of runNumbers
% startDigitsArray = [];  % Array of startDigits (or empty array for individual digits)
% endDigitsArray = [];  % Array of endDigits (or empty array for individual digits)
individualDigitsArray = {
    [001, 006, 007, 008, 011, 012, 014, 016, 019, 020], 
    [003, 004, 005, 006, 007, 008, 009, 010, 011, 012], 
    [001, 007, 010, 013, 016, 019, 021, 022, 023, 024], 
    [001, 004, 005, 007, 008, 010, 011, 012, 014, 016]};  % Cell array of individual digits for each set

[medians,  percentiles25, percentiles75] = calculateTemporalStats(runNumbers, individualDigitsArray);

% Multiply all values in medians and percentiles75 by 1000
for setIndex = 1:numel(runNumbers)
    medians{setIndex} = medians{setIndex} * 1000;
    percentiles25{setIndex} = percentiles25{setIndex} * 1000;
    percentiles75{setIndex} = percentiles75{setIndex} * 1000;
end

% plotting
customColors = [
    0 0 0;        % Black
    0.5 0.5 0.5;  % Light Gray
    1 0 0;        % Red
    1 0.6 0.6;    % Light Red
];

% Define x-axis tick values
xTicks = [0, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000];

% Plotting
figure;
hold on;

for setIndex = 1:numel(runNumbers)
    colorIndex = mod(setIndex - 1, size(customColors, 1)) + 1;
    setColor = customColors(colorIndex, :);
    xValues = xTicks(1:numel(medians{setIndex}));

    % Plot median data with thicker lines, markers, and legend
    medianPlot = plot(xValues, medians{setIndex}, 'o-', 'Color', setColor, 'LineWidth', 2, 'MarkerFaceColor', setColor, 'DisplayName', sprintf('Set %d Median', setIndex));

    % Plot the area between 25th and 75th percentiles with matching fill and dashed outline
    fill([xValues, fliplr(xValues)], [percentiles25{setIndex}, fliplr(percentiles75{setIndex})], setColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    plot(xValues, percentiles25{setIndex}, 'Color', setColor, 'LineStyle', '--', 'HandleVisibility', 'off');
    plot(xValues, percentiles75{setIndex}, 'Color', setColor, 'LineStyle', '--', 'HandleVisibility', 'off');

    % Add legend entry for median data only
    legendEntries{setIndex} = medianPlot;
end

% Add legend
legend([legendEntries{1}, legendEntries{2}, legendEntries{3}, legendEntries{4}], ...
  'Run 1, Q = 0.25 L/s, S = 0.01', 'Run 2, Q = 0.12 L/s, S = 0.01', 'Run 3, Q = 0.25 L/s, S = 0.02', 'Run 4, Q = 0.12 L/s, S = 0.02');

% Add labels and title
xlabel('Dimensionless time, {\it t*}');
ylabel('Channel depth, {\it h_{channel}} (mm)')

% Set x-axis tick values
xticks(xticks);
xticklabels(arrayfun(@(x) sprintf('%g', x), xticks/1000, 'UniformOutput', false));

% Remove top border
ax = gca; % Get the current axes
ax.Box = 'on'; % Turn off the box around the axes

%% calculate channel aggradation 
runNumbers = [1, 2, 3, 4];  % Array of runNumbers

individualDigitsArray = {
    [007, 032], 
    [008, 013], 
    [006, 025], 
    [006, 020]};  % Cell array of individual digits for each set

medians = differenceDEM(runNumbers, individualDigitsArray);
timeHours = zeros(size(individualDigitsArray));
for i = 1:numel(individualDigitsArray)
    difference = individualDigitsArray{i}(2) - individualDigitsArray{i}(1);
    timeHours(i) = difference;
end
aggrad_rate=(medians(:)./timeHours).*1000;

%normalized to flow depth
flow_depth_averages=[0.01;0.0075;0.004;0.0031].*1000;
norm_aggrad_rate=aggrad_rate./flow_depth_averages


%% Figure 7: plot empirical pdfs of scour depth and flow depth

M=[M1 M2 M3 M4 M5].';
%fit scour depths to gamma
Pd = fitdist(M,'Gamma');
[heights, locations] = hist(M, 25);
width=locations(2)-locations(1);
heights=heights/(length(M)*width); 
l=linspace(0, max(M), 1000);
P=pdf(Pd,l);

% plot scour depths for each run 
figure(1)
subplot(1,2,2)
hold on
%plot(l.*1000, P, 'LineWidth', 2.5, 'color', 'black')
%plot(l.*1000, P, 'LineWidth', 2.5, 'color', [0.5 0.5 0.5])
%plot(l.*1000, P, 'LineWidth', 2.5, 'color', [1 0 0])
plot(l.*1000, P, 'LineWidth', 2.5, 'color', [1 0.6 0.6])

xlabel('Modeled flow depths (mm)') 
ylabel('Probability density')
legend('Run 1, Q = 0.25 L/s, S = 0.01', 'Run 2, Q = 0.12 L/s, S = 0.01', 'Run 3, Q = 0.25 L/s, S = 0.02', 'Run 4, Q = 0.12 L/s, S = 0.02')
%% Figure 9 - fitting set thicknesses to exponential
% calculate Paola & Borgman pdf
% beta = 1.9;
% l = linspace(0, 25, 16);
% a = 1/beta; 
% f = @(s)((a.*exp((-1).*a.*s)).*(exp((-1).*a.*s)+(a.*s)-1))./((1-exp((-1).*a.*s)).^2); % from Paola and Borgman
% setPB=f(l); 
% setPB(1)=0; % first value comes out as NaN

% calculate Paola & Borgman pdf
mean_set=2.8; % observed
bins=26;
a_mean = 1.645/mean_set;
l = linspace(0, 25, bins);
f = @(s)((a_mean.*exp((-1).*a_mean.*s)).*(exp((-1).*a_mean.*s)+(a_mean.*s)-1))./((1-exp((-1).*a_mean.*s)).^2); % from Paola and Borgman
setPB=f(l); % pdf for predicted 
setPB(1)=0; % first value comes out as NaN

% histogram of observations
binCenters = [l(1:end) - 0.5 * l(2), l(end)]; % add another value at the very end
set_obs=histcounts(setRun4.*1000, binCenters, 'Normalization', 'probability');

%plotting
figure
bar(l, set_obs, 'FaceColor', 'white');
hold on
plot(l, setPB/sum(setPB), 'LineWidth', 2, 'Color', 'blue');
legend('Observed sets', 'Predicted PDF')
xlabel('set thickness,{\it s} (mm)')
ylabel('Probability density')
title('Run 4, Q = 0.12 L/s, S = 0.02');


% goodness-of-fit test
chi_2=sum((set_obs(2:end)-setPB(2:end)).^2/setPB(2:end));

% Degrees of freedom
df = numel(l) - 1;

% P-value from chi-squared distribution
pValue = 1 - chi2cdf(chi_2, df);

% Display results
fprintf('Chi-squared test results:\n');
fprintf('Chi2 statistic = %.4f\n', chi_2);
fprintf('Degrees of freedom = %d\n', df);
fprintf('P-value = %.4f\n', pValue);

% Determine whether to reject or accept the null hypothesis
alpha = 0.05; % Set your significance level
if pValue < alpha
    fprintf('Reject the null hypothesis. The distributions are different.\n');
else
    fprintf('Accept the null hypothesis. The distributions are similar.\n');
end

cv_sets=std(setRun4)/mean(setRun4);

%% Figure 10 
% A: set thickness violin plots
% note: variables generated from this part of code are saved as
% violinplt_scour.mat and violinplt_set.mat
pd1 = fitdist(Zclean,'Kernel');
y1 = random(pd1,10000,1);
y1(end+1)=max(Zclean);
y1(end+1)=min(Zclean);

s=size(y4,1)
origin4=cell(s,1);
origin4(1:end)={'Run 4'};

% run from here to replot
set=[y1;y2;y3;y4]
origin=[origin1;origin2;origin3;origin4]
grouporder={'Run 1','Run 2','Run 3','Run 4'};

% plot
figure
violins=violinplot(set.*1000, origin, 'GroupOrder', grouporder, 'ShowData', false, 'ViolinAlpha', 0.6,'MedianColor', [0 0 0], 'BoxColor', [0 0 0],'EdgeColor', [0 0 0],'ShowNotches', true)
violins(1).ViolinColor = [0.2 0.2 0.2]; % run 1 dark gray
violins(2).ViolinColor = [0.8 0.8 0.8]; % run 2 gray 
violins(3).ViolinColor = [1 0 0]; % run 3 red 
violins(4).ViolinColor = [1 0.6 0.6] % run 4 pink
ylabel('set thickness (mm)')

% B: plot and calculate preservation ratio
customColors = [
    0 0 0;        % Black
    0.5 0.5 0.5;  % Light Gray
    1 0 0;        % Red
    1 0.6 0.6;    % Light Red
];
x_axes=["Run 1","Run 2","Run 3","Run 4"];
C = categorical(x_axes);
PSRatios=[0.3824, 0.3135, 0.4939, 0.6173];
figure
hold on
for i = 1:length(C)
    scatter(C(i), PSRatios(i), 100, customColors(i, :), 's', 'MarkerFaceColor', customColors(i, :))
end

ylabel('preservation ratio')
legend('Q = 0.25 L/s, S = 0.01', 'Q = 0.12 L/s, S = 0.01', 'Q = 0.25 L/s, S = 0.02', 'Q = 0.12 L/s, S = 0.02')
ax = gca; % Get the current axes
ax.Box = 'on'; % Turn off the box around the axes
ylim([0.2 0.7])


%% calculating set thickness for braiding 
setRun = calculateSet(1, [007, 032]) % Run 1

%% plotting the strat record
%Function to load Robert's data and visualize 
%load data

% converting to Roberts data demo 
timeHours=cell2mat(timeHours');
horz_raw=cell2mat(horz_raw);
Depth=horz_raw.';
xPosclip=loc_y{1}; 
time = timeHours(2:end);
x = abs(xPosclip-xPosclip(1));
Z = Depth;
Z = fliplr(Z); %comment out for working with Guala data
dx = mean(diff(x));
dt = mean(diff(time));
nt=numel(time);
newRow=-1.555*ones(1, size(Z,2)); % add new row 
Z=vertcat(newRow,Z);

% plot section
nt=numel(time);
smpInterval = 1;
figure('Position',[100,100,800,1200]);
subplot(4,1,4)
[S,Z_braid]=plot_Section(Z,x,smpInterval);
hold on 
plot(x, Z_braid(5,:),'LineWidth', 2, 'Color', 'k')
title('Run 4')
%xlabel('Cross Stream Distance (m)')
%ylabel('Stratigraphic Surface Elevation (m)')
xlim([0.7 1.8])
ylim([-1.425 -1.39])
mycolormap = customcolormap([0 .25 .5 .75 1], {'#9d0142','#f66e45','#ffffbb','#65c0ae','#5e4f9f'});
cb=colorbar();
colormap(mycolormap);
cb.YTick=[1 size(Zlid,2)];
dimt=round(((15*3600)*9.8*(0.00042^2))/(0.00012*(0.02^(2/3))));
cb.YTickLabel={'0', num2str(dimt)}; %dimensionless time 
title(cb, 't*')
%saveas(gcf,'setsections.fig')
[ax1,h1]=suplabel('Cross Stream Distance (m)');
[ax2,h2]=suplabel('Stratigraphic Surface Elevation (m)','y');
% plot section 
figure
subplot(4,1,1)
[F,hsrf,hlyr]=section(S(1:smpInterval:end,:),[],x);
subplot(2,1,2)
[F,hsrf]=sectioncol(S(1:smpInterval:end,:),[],x);

x1=channelbelt1(end);
x2=channelbelt2(end);
% plot set thickness  
diff_S=S(:,x1:x2);
diff_S=diff(diff_S);
diff_S(diff_S == 0) = NaN;
figure
boxplot(diff_S.')
set(gca, 'fontsize', 15)
ylabel('set thickness (m)')
xlabel('time (hours)') 
%% plot top 75 percent of scours

[file_list, path_n]= uigetfile('.mat', 'grab the files', 'MultiSelect', 'on')

% get the file paths
if iscell(file_list) == 0 
    file_list = (file_list); 
end

data_in=load([path_n, file_list], 'data');
loc_x=data_in.data.xVec; 
loc_y=data_in.data.yVec;
loc_z=data_in.data.elevationMasked;
loc_zd=data_in.data.elevationMaskedDetrended;
loc_zraw=data_in.data.elevationRaw;
clean=find(all(isnan(loc_z),1)); %columns with all NaNs
loc_zraw(:,clean)=[];
loc_zd(:,clean)=[];
loc_z(:,clean)=[];
loc_x(clean)=[];

[Zlid] = scourspace(loc_z,loc_zraw,loc_x);
Zlid(isnan(Zlid)) = [];
Zlid75(isnan(Zlid75))=0;
Zlid75(Zlid75<0)=1;
blue = cat(3, zeros(size(loc_zd)), zeros(size(loc_zd)), zeros(size(loc_zd))+0.5);

figure(1)
a(1)=subplot(4,2,5)
b=imagesc(loc_x,loc_y,loc_zd.*1000)
pos1 = get(a,'Position');
colormap gray
caxis([-5 5])
set(b, 'AlphaData', ~isnan(loc_zd))
hold on 
h=imagesc(loc_x, loc_y, blue)
set(h, 'AlphaData', Zlid75)
hold off
ylabel('Distance (m)')
xlim([18 34])
ylim([0.5 2.5])
%c=colorbar('Northoutside')
%c.Label.String='Detrended Elevation (mm)'
set(a(1),'Position',pos1)
%xlabel('Distance (m)')

subplot(4,2,2)
b=imagesc(loc_x,loc_y,loc_zd.*1000)
colormap gray
caxis([-5 5])
set(b, 'AlphaData', ~isnan(loc_zd))
hold on 
h=imagesc(loc_x, loc_y, blue)
set(h, 'AlphaData', Zlid75)
hold off
ylabel('Distance (m)')
xlim([18 34])
ylim([0.5 2.5])

subplot(4,1,3)
b=imagesc(loc_x,loc_y,loc_zd.*1000)
colormap gray
caxis([-5 5])
set(b, 'AlphaData', ~isnan(loc_zd))
hold on 
h=imagesc(loc_x, loc_y, blue)
set(h, 'AlphaData', Zlid75)
hold off
ylabel('Distance (m)')
xlim([18 34])
ylim([0.5 2.5])

subplot(4,1,4)
b=imagesc(loc_x,loc_y,loc_zd.*1000)
colormap gray
caxis([-5 5])
set(b, 'AlphaData', ~isnan(loc_zd))
hold on 
h=imagesc(loc_x, loc_y, blue)
set(h, 'AlphaData', Zlid75)
hold off
ylabel('Distance (m)')
xlim([18 34])
ylim([0.5 2.5])
xlabel('Distance (m)')