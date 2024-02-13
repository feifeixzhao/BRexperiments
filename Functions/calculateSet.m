function setRun = calculateSet(run_number, time_range)
    folder_name = ['run', num2str(run_number)];
    mat_files = dir(fullfile(folder_name, '*.mat'));
    % Filter mat_files based on time_range
    selected_files = {};
    for i = 1:length(mat_files)
        filename = mat_files(i).name;
        time_match = regexp(filename, 'processedElevation_(\d+)', 'tokens');
        if ~isempty(time_match)
            file_time = str2double(time_match{1}{1});
            if file_time >= time_range(1) && file_time <= time_range(2)
                selected_files{end+1} = fullfile(folder_name, filename); % Include full path
            end
        end
    end
    num_files = length(selected_files);
    display(selected_files)

        % Preallocate cell arrays
        data_in = cell(1, num_files);
        loc_x = cell(1, num_files);
        loc_y = cell(1, num_files);
        loc_z = cell(1, num_files);
        loc_zd = cell(1, num_files);
        loc_zraw = cell(1, num_files);
        timeHours = zeros(1, num_files);
        clean = cell(1, num_files);
        channelbelt1 = cell(1, num_files);
        channelbelt2 = cell(1, num_files);
        horz = cell(1, num_files);
        horz_raw = cell(1, num_files);
        Z = cell(1, num_files);
        S = cell(1, num_files);
        diff_S = cell(1, num_files);
        sthickness = cell(1, num_files);
        
        for i = 1:num_files
            filename = selected_files{i};
            filename = char(filename);
            data_in{i} = load(selected_files{i}, 'data');
            loc_x{i} = data_in{i}.data.xVec; 
            loc_y{i} = data_in{i}.data.yVec;
            loc_z{i} = data_in{i}.data.elevationMasked;
            loc_zd{i} = data_in{i}.data.elevationMaskedDetrended;
            loc_zraw{i} = data_in{i}.data.elevationRaw;  
            timeHours(i) = data_in{i}.data.metadata.timeHours;
            clean{i} = all(isnan(loc_z{i}), 1);
            loc_z{i}(:, clean{i}) = [];
        end
        
        for j = 1:num_files
            loc_zraw{j}(:, clean{j}) = [];
            loc_x{j}(clean{j}) = [];
        end
        
        % preparing for Roberts code
        time = timeHours(2:end);
        xPosclip = loc_y{1}; 
        x = abs(xPosclip - xPosclip(1));
        dx = mean(diff(x));
        dt = mean(diff(time));
        nt = numel(time);
        
        for m = 1:length(loc_z{1})
            [channelbelt1{m}, channelbelt2{m}, horz{m}, horz_raw{m}] = scourtime(m, loc_z, loc_zraw);
        end
        
        horz = cellfun(@cell2mat, horz, 'UniformOutput', false);
        horz_raw = cellfun(@cell2mat, horz_raw, 'UniformOutput', false);
        Z = cellfun(@(a) a.', horz_raw, 'UniformOutput', false);
        Z = cellfun(@fliplr, Z, 'UniformOutput', false); 
        
        % calculate strat from topo  
        for section_idx = 1:size(Z,2)
                % Calculate S for each cell in Z
                calculated_S = plot_Section(Z{section_idx});
                S{section_idx}=calculated_S;
        end
        
        diff_S=cellfun(@diff, S, 'UniformOutput', false); 
        
        % Remove the first cell from loc_z
        loc_z(1) = [];
        
        % Convert diff_S into a 1x2 cell with 1351x17501 arrays
        diff_S_new = cell(1, size(loc_z, 2));
        for i = 1:size(loc_z, 2)
            diff_S_new{i} = reshape(cell2mat(cellfun(@(x) x(i, :), diff_S, 'UniformOutput', false)), 1351, 17501);
        end

    
    % Mask out NaN values in diff_S_new using loc_z and set zeros to NaN
    masked_diff_S = cell(size(diff_S_new));
    for i = 1:numel(diff_S_new)
        mask = ~isnan(loc_z{i});
        masked_diff_S{i} = diff_S_new{i};
        masked_diff_S{i}(~mask | masked_diff_S{i} <= 0.001) = NaN; % removing below 1mm
    end
    
    setRun = cell2mat(masked_diff_S);
end


