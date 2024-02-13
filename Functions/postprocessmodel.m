function M=postprocessmodel(run_number, depth_time_number, elevation_time_number)

    % Specify the base directory for depths and elevations
    base_depth_data_dir = 'C:\Users\feife\Box\BR_experiments\caesar_lisflood_2022\new_flow_depths\';
    base_elevation_data_dir = 'C:\Users\feife\Box\BR_experiments\braidingExperimentsTopo_v2\';
    
    % Specify the run number
    run_str = sprintf('run%d', run_number); % Construct the run folder name
    
    % Specify the time number for depths
    depth_time_str = sprintf('%03d', depth_time_number); % Format depth time number as two digits
    
    % Construct the full path to depth data files
    depthfile1 = fullfile(base_depth_data_dir, run_str, [run_str, '_', depth_time_str, 'hr_waterdepth60', '.txt']);
    depthfile2 = fullfile(base_depth_data_dir, run_str, [run_str, '_', depth_time_str, 'hr_waterdepth60','.asc']);
    depthfile3 = fullfile(base_depth_data_dir, run_str, [run_str, '_', depth_time_str, 'hr_waterdepth60','.mat']);

    load(depthfile3, 'M')
    
    % Specify the time number for elevations
    elevation_time_str = sprintf('%03d', elevation_time_number); % Format elevation time number as three digits
    
    % Construct the full path to the processedElevation .mat file
    elevation_mat_path = fullfile(base_elevation_data_dir, run_str, [run_str, '_processedElevation_', elevation_time_str, '.mat']);
    load(elevation_mat_path, 'data'); % Load elevation data from the constructed .mat file path

    %% script to mask out the values outside the channel belt

    options.dx = 0.020677; % original DEM grid spacing
    options.crop_lims.x = [0 36.42]; % crop original DEM using x1 x2 y1 y2, in meters. Do this, rather than row-column, to account for different DEM extents
    options.crop_lims.y = [0.1 2.6];
    options.downsample_factor = 10; % i.e., downsample by this factor
    options.add_sidewalls = true;
    downsample_factor = options.downsample_factor; % i.e., downsample input DEM by this factor

    [dem.x, dem.y] = meshgrid(data.xVec, data.yVec);
    dem.z = data.elevationMasked; % or elevationMasked

    % catch the NoData value and reassign to NaN
    dem.z(dem.z == -9999) = NaN;

    % check that the elevations and coordinates are in m, not mm, using the range in elevations
    if range(dem.z(:)) > 1
        dem.z = dem.z / 1000;
        dem.x = dem.x / 1000;
        dem.y = dem.y / 1000;
    end

    rows = (1:size(dem.z, 1))';
    cols = 1:size(dem.z, 2);

    % crop using input (x, y) coordinates by finding the corresponding
    % (column, row) coordinates in DEM
    temp = interp1(dem.x(1, :), cols, options.crop_lims.x, 'nearest');
    options.crop_lims.col_start = temp(1);
    options.crop_lims.col_end = temp(2);
    temp = interp1(dem.y(:, 1), rows, options.crop_lims.y, 'nearest');
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
    if size(dem.z, 1) < size(dem.z, 2)
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
    if downsample_factor > 1
        dem.z = imresize(dem.z, 1 / downsample_factor);
        dem.x = imresize(dem.x, 1 / downsample_factor);
        dem.y = imresize(dem.y, 1 / downsample_factor);
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
    Mask.grid = dem.z - min(dem.z(:)); % set lowest elevation as zero
    Mask.x = dem.x(1, :);
    Mask.y = dem.y(:, 1);
    dx = abs(diff(Mask.x));
    dy = abs(diff(Mask.y));
    Mask.dx = dx(1);
    Mask.dy = dy(1);
    clear dem

    % put mask on DEM
    M.grid(isnan(Mask.grid)) = nan;
    M.grid(M.grid<0.001) = nan;
end

