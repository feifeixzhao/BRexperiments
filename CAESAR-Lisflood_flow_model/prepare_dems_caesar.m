function prepare_dems_caesar(infile,outfile,options)
% This function modifies DEMs for use in the CAESAR-LISFLOOD flow model.
% Created January 27, 2017 by Ajay Limaye (aslimaye@umn.edu).
% Last edited February 8, 2017 by Ajay Limaye (aslimaye@umn.edu).

downsample_factor = options.downsample_factor; % i.e., downsample input DEM by this factor

load(infile,'data')
[dem.x,dem.y] = meshgrid(data.xVec,data.yVec);
dem.z = data.elevationRaw; % or elevationMasked

% catch the NoData value and reassign to NaN
dem.z(dem.z==-9999)=NaN;

% % check that the elevations and coordinates are in m, not mm, using the range in elevations
% if range(dem.z(:))>1
%     dem.z = dem.z/1000;
%     dem.x = dem.x/1000;
%     dem.y = dem.y/1000;
% end

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

% Fill any holes in dem using the first non-NaN elevation in
% the row (CAESAR can't accept NoData values in DEM)
[nan_i,nan_j] = find(isnan(dem.z));
for i=1:numel(nan_i)
    nan_j2 = find(~isnan(dem.z(nan_i(i),:)),1,'first');
    dem.z(nan_i(i),nan_j(i)) = dem.z(nan_i(i),nan_j2);
end

% resample
if downsample_factor>1
dem.z = imresize(dem.z,1/downsample_factor);
dem.x = imresize(dem.x,1/downsample_factor);
dem.y = imresize(dem.y,1/downsample_factor);
end

% add a shallow moat to north side
% moat_depth = 0.1*std(dem.z(1,:)); % if too deep, takes too long to fill; if too shallow, doesn't spread in moat enough before entering DEM
% moat_width = ceil(0.01*size(dem.z,1)); % arbitrary
% wall_height = 1; % 1m, i.e., high enough nothing will spill out
% sidewall_length = moat_width + 1;

%%% This is the setup for a moat with no spillover
%z_add = [repmat(max(dem.z(:))+wall_height,1,size(dem.z,2));repmat(min(dem.z(1,:))-moat_depth,moat_width,size(dem.z,2))];

% %%% This is the setup for a spillover
% z_add = [repmat(max(dem.z(:))+wall_height,1,size(dem.z,2));...
%                 repmat(max(dem.z(1,:))-moat_depth,moat_width,size(dem.z,2));...
%                 repmat(max(dem.z(1,:)),1,size(dem.z,2))]; % this puts a lip in place so water can get to all of dem
% 
% % add sidewalls for moat section
% z_add(1:sidewall_length,1)=max(dem.z(:))+wall_height;
% z_add(1:sidewall_length,end)=max(dem.z(:))+wall_height;
% dem.z = [z_add;dem.z];

% if specifed, add sidewalls the span the full domain 
if options.add_sidewalls
    dem.z(:,[1 end])=max(dem.z(:));
end

% check for NoData values
if any(or(isnan(dem.z(:)),dem.z(:)<-1000))
    error('DEM has NoData values')
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
M.grid = dem.z-min(dem.z(:)); % set lowest elevation as zero
M.x = dem.x(1,:); 
M.y = dem.y(:,1);
dx=abs(diff(M.x));
dy=abs(diff(M.y));
M.dx=dx(1);
M.dy=dy(1);
clear dem

imagesc(M.x(:),M.y(:),M.grid)
xlabel('Distance (m)')
ylabel('Distance(m)')
cbar = colorbar;
set(get(cbar,'title'),'string','Elevation (m)')

WriteArcGrid(M,[outfile,'.asc'])
% geotiffwrite(outfile, dem.z, R_dem, ...
%             'GeoKeyDirectoryTag', info_dem.GeoTIFFTags.GeoKeyDirectoryTag);
movefile([outfile,'.asc'],[outfile,'.txt']) % CAESAR requires a .txt extension
sprintf('Wrote DEM to %s',[outfile,'.txt'])
end