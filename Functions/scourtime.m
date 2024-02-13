function [channelbelt1, channelbelt2, horz, horz_raw] = scourtime(interval, Z, Z_raw)
    horz = cellfun(@(a) a(:, interval), Z, 'UniformOutput', false);
    horz_raw = cellfun(@(a) a(:, interval), Z_raw, 'UniformOutput', false);
    channelbelt1 = nan(1, length(horz));
    channelbelt2 = nan(1, length(horz));
    
    for i = 1:length(horz) % plotting profiles 
        non_nan_indices = find(~isnan(horz{i}));
        
        if ~isempty(non_nan_indices)
            channelbelt1(i) = non_nan_indices(1);
            channelbelt2(i) = non_nan_indices(end);
        end
    end
end
