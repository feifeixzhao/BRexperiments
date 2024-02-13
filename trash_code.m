%% compute beta values spatially - all the times combined
set_onem=reshape(diff_S, 37, []); % for every one meter 
set_3m=reshape(diff_S, 1591, []); % for ever 3.4 meters
% calculate beta every 1 meter 
newA=cellfun(@(c) c(c>=0.001), set_onem, 'UniformOutput', false);
newA2=arrayfun(@(x) vertcat( newA{:,x}), 1:size(newA,2), 'Uniform', 0) %resized 
% calculate beta every 3.4 meter 
new3m=cellfun(@(c) c(c>=0.001), set_3m, 'UniformOutput', false);
new3m2=arrayfun(@(x) vertcat(new3m{:,x}), 1:size(new3m,2), 'Uniform', 0) %resized 

avg_newA=cell(size(newA2)); 
avg_new3m=cell(size(new3m2));
for i=1:numel(newA2)
    avg_newA{i}=mean(newA2{i}); %mean a value 
    std_newA{i}=std(newA2{i}); %std a value
end
for i=1:numel(new3m2)
    avg_new3m{i}=mean(new3m2{i});
    violin{i}=fitdist(new3m2{i},'Kernel');
    yviolin{i}=random(violin{i},10000,1);
    yviolin{i}(end+1)=max(new3m2{i});
    yviolin{i}(end+1)=min(new3m2{i});
end

new=cell2mat(yviolin); % for violin plots

beta_1m=cell2mat(avg_newA)./1.64493; % formula for beta at 1 meter
a_beta_1m=1./beta_1m;
betastd_1m=cell2mat(std_newA)./1.64493; % std value for beta at 1 meter 
beta_3m=cell2mat(avg_new3m)./1.64493;

%%
y = num2cell(1:numel(Zclean));
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], Zclean, y, 'UniformOutput', 0); % adding labels to the cells
X = vertcat(x{:});
boxplot(X(:,1), X(:,2), 'labels',{'t = 5 hrs','t = 15 hrs', 't = 25 hrs', 't = 30 hrs', 't = 35 hrs'})




function [Zlid] = scourspace1(Z, Z_raw, x)
% Calculates the cross-stream scour profile across space of Ajay's braided river DEMs.
% Written by Feifei Zhao, 2021.

% Preallocate variables
channelbelt1 = nan(1, size(Z, 2));
channelbelt2 = channelbelt1;
data1 = nan(51, size(Z, 2));
data2 = data1;

% Find channel belts and extract data
for i = 1:size(Z, 2)
    index1 = find(~isnan(Z(:, i)), 1, 'first');
    index2 = find(~isnan(Z(:, i)), 1, 'last');

    if ~isempty(index1)
        channelbelt1(i) = index1; 
        channelbelt2(i) = index2;
    end
    data1(:, i) = Z_raw(channelbelt1(i)-50:channelbelt1(i), i);
    data2(:, i) = Z_raw(channelbelt2(i):channelbelt2(i)+50, i);
end

% Compute lid and interpolate
data = [data1; data2];
Lid=mean(data, 1); % lid now has multiple values (aka slanted lid)
% interpolate
interp = linspace(Lid(1), Lid(end), length(Lid));

% subtracting the lid
Zlid=interp-Z;
end

%% scour depths and flow depths 

%combine data with scour depths
data=[y4;MaskedM_all];

%index data
s1=size(y4,1);
origin1=cell(s1,1);
origin1(1:end)={'Set Thickness'};
s2=size(MaskedM_all,1);
origin2=cell(s2,1);
origin2(1:end)={'Flow Depths'};

origin=[origin1;origin2]
grouporder={'Set Thickness','Flow Depths'};

figure 
subplot(1,4,1)
violins=violinplot(data.*1000,origin, 'GroupOrder', grouporder, 'ShowData', false, 'ViolinAlpha', 0.6,'MedianColor', [0 0 0], 'BoxColor', [0 0 0],'EdgeColor', [0 0 0],'ShowNotches', true);
violins(1).ViolinColor = [0.2 0.2 0.2]; % run 1 dark gray
violins(2).ViolinColor = [1 0 0]; % red
ylim([-1 50])
ylabel('mm')

subplot(1,4,2)
violins=violinplot(data.*1000,origin, 'GroupOrder', grouporder, 'ShowData', false, 'ViolinAlpha', 0.6,'MedianColor', [0 0 0], 'BoxColor', [0 0 0],'EdgeColor', [0 0 0],'ShowNotches', true);
violins(1).ViolinColor = [0.2 0.2 0.2]; % run 1 dark gray
violins(2).ViolinColor = [1 0 0]; % red
ylim([-1 50])
ylabel('mm')

subplot(1,4,3)
violins=violinplot(data.*1000,origin, 'GroupOrder', grouporder, 'ShowData', false, 'ViolinAlpha', 0.6,'MedianColor', [0 0 0], 'BoxColor', [0 0 0],'EdgeColor', [0 0 0],'ShowNotches', true);
violins(1).ViolinColor = [0.2 0.2 0.2]; % run 1 dark gray
violins(2).ViolinColor = [1 0 0]; % red
ylim([-1 50])
ylabel('mm')

subplot(1,4,4)
violins=violinplot(data.*1000,origin, 'GroupOrder', grouporder, 'ShowData', false, 'ViolinAlpha', 0.6,'MedianColor', [0 0 0], 'BoxColor', [0 0 0],'EdgeColor', [0 0 0],'ShowNotches', true);
violins(1).ViolinColor = [0.2 0.2 0.2]; % run 1 dark gray
violins(2).ViolinColor = [1 0 0]; % red
ylim([-1 50])
ylabel('mm')

function [channelbelt1,channelbelt2,horz,horz_raw,lid,Zlid] = scourtime(interval,Z,Z_raw)
% calculates the cross stream scour profile across different hours of Ajay's braided river DEMs
% written by Feifei Zhao, 2021
horz=cellfun(@(a) a(:,interval), Z, 'UniformOutput', false);
horz_raw=cellfun(@(a) a(:,interval), Z_raw, 'UniformOutput', false);
channelbelt1=nan(1, size(Z, 2));
channelbelt2=channelbelt1;

    for i=1:length(horz) %plotting profiles 
        channelbelt1(i) = find(~isnan(horz{i}), 1, 'first');
        channelbelt2(i) = find(~isnan(horz{i}), 1, 'last');
        lid(i)=floodplain1(horz_raw{i}, channelbelt1(i), channelbelt2(i));
    end
    
lid=mean(lid);
Zlid=cellfun(@(x) lid-x, horz, 'un', 0);
end


% Script to convert CAESAR output files for water depth to .mat format

clear,clc
addpath('source')
data_dir = 'C:\Users\feife\Box\BR_experiments\caesar_lisflood_2022\new_flow_depths\run4\';
times = 60;

for i=1:numel(times)
    depthfile1 = [data_dir,'run4_011hr_waterdepth',num2str(times(i)),'.txt'];
    depthfile2 = [data_dir,'run4_011hr_waterdepth',num2str(times(i)),'.asc'];
    depthfile3 = [data_dir,'run4_011hr_waterdepth',num2str(times(i)),'.mat'];
    if and(exist(depthfile1),~exist(depthfile3))
       copyfile(depthfile1,depthfile2)
        M = ReadArcGrid(depthfile2);
        save(depthfile3,'M')
    end
end

load(depthfile3,'M')

%% script to mask out the values outside channel belt

options.dx = 0.020677; % original DEM grid spacing
options.crop_lims.x = [0 36.42]; % crop original DEM using x1 x2 y1 y2, in meters. Do this, rather than row-column, to account for different DEM extents
options.crop_lims.y = [0.1 2.6];
options.downsample_factor = 10; % i.e., downsample by this factor
options.add_sidewalls = true;
downsample_factor = options.downsample_factor; % i.e., downsample input DEM by this factor


load('run4_processedElevation_012.mat');
[dem.x,dem.y] = meshgrid(data.xVec,data.yVec);
dem.z = data.elevationMasked; % or elevationMasked

% catch the NoData value and reassign to NaN
dem.z(dem.z==-9999)=NaN;

% check that the elevations and coordinates are in m, not mm, using the range in elevations
if range(dem.z(:))>1
    dem.z = dem.z/1000;
    dem.x = dem.x/1000;
    dem.y = dem.y/1000;
end

rows = (1:size(dem.z,1))';
cols = 1:size(dem.z,2);

% crop using input (x,y) coordinates by finding the corresponding
% (column,row) coordinates in DEM
temp = interp1(dem.x(1,:),cols,options.crop_lims.x,'nearest');
options.crop_lims.col_start = temp(1);
options.crop_lims.col_end = temp(2);
temp = interp1(dem.y(:,1),rows,options.crop_lims.y,'nearest');
options.crop_lims.row_start = temp(1);
options.crop_lims.row_end = temp(2);

dem.z = dem.z(options.crop_lims.row_start:options.crop_lims.row_end, options.crop_lims.col_start:options.crop_lims.col_end); % crop dem
dem.x = dem.x(options.crop_lims.row_start:options.crop_lims.row_end, options.crop_lims.col_start:options.crop_lims.col_end); % crop dem
dem.y = dem.y(options.crop_lims.row_start:options.crop_lims.row_end, options.crop_lims.col_start:options.crop_lims.col_end); % crop dem

% Get the size of your DEM
[rows, cols] = size(dem.z);

% Define the default value
default_value = 1;

% Loop through each column
for col = 1:cols
    % Find the NaN values in the current column
    nan_indices = isnan(dem.z(:, col));
    
    if any(~nan_indices)
        % Find the first and last non-NaN values in the column
        first_non_nan_index = find(~nan_indices, 1, 'first');
        last_non_nan_index = find(~nan_indices, 1, 'last');
        
        % Replace NaN values between the first and last non-NaN values with the default value
        dem.z(first_non_nan_index+1:last_non_nan_index-1, col) = default_value;
    end
end

% Orient the DEM so that the long direction is north-south, which the
% following steps assume. The DEM gets rotated back to its original
% orientation at the end.
rotated_dem = false;
if size(dem.z,1)<size(dem.z,2)
    rotated_dem = true;
    dem.z = dem.z';
    % the x coordinates and y coordinates switch
    tempx = dem.x;
    tempy = dem.y;
    dem.x = tempy;
    dem.y = tempx;
    clear tempx tempy
end


% resample
if downsample_factor>1
dem.z = imresize(dem.z,1/downsample_factor);
dem.x = imresize(dem.x,1/downsample_factor);
dem.y = imresize(dem.y,1/downsample_factor);
end

% rotate back to original coordinates
if rotated_dem
    dem.z = dem.z';
    % the x coordinates and y coordinates switch
    tempx = dem.y;
    tempy = dem.x;
    dem.x = tempx;
    dem.y = tempy;
    clear tempx tempy
end

% format as a structure array for WriteArcGrid
Mask.grid = dem.z-min(dem.z(:)); % set lowest elevation as zero
Mask.x = dem.x(1,:); 
Mask.y = dem.y(:,1);
dx=abs(diff(Mask.x));
dy=abs(diff(Mask.y));
Mask.dx=dx(1);
Mask.dy=dy(1);
clear dem

% put mask on DEM
M.grid(isnan(Mask.grid)) = nan;
M.grid(M.grid<0.001) = nan

%%
% Load your data
load('/Users/feife/Box/BR_experiments/braidingExperimentsTopo_v2/run1/run1_processedElevation_013.mat'); % Make sure to specify the correct file path

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



