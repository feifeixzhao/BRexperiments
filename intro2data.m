% readme.m: Readme file for time series data for experimental braided channels, including
% example plotting scripts. Based on data from Limaye (2020), How do
% braided rivers grow channel belts?, JGR-Earth Surface,
% doi:10.1029/2020JF005570.

%%% Overview 
% The data are organized in separate directories for each of four runs
% of the experiment. Each file name includes the name of the run
% (run1,run2,...) and a number to indicate the sequential order of the 
% topography measurements (001,002,...).

%%% Variables
% In each file, 'data' is a MATLAB structure array with the following
% fields:
% - xVec: vector with x-coordinates of elevation data. x = 0 is the inlet
% of the basin. 
% - yVec: vector with y-coordinates of elevation data. 
% - elevationMasked: DEM, with elevations masked by the extent of the
% channel belt. Areas outside the channel belt are set to NaN.
% - elevationMaskedDetrended: DEM, with elevations masked by the extent of
% the channel belt and after a linear trend has been removed. 
% - metadata: Includes metadata for this observation stored in separate
% structure fields:
% -- timeHours: elapsed experiment run time, in hours
% -- Q: water discharge (m^3/s)
% -- S: bed slope (dimensionless)
% -- runName: name of the experiment run. 

%%% Plotting examples

%% to plot *detrended* elevation masked by the extent of the channel belt
figure
tiledlayout(3,1)
nexttile
imagesc(data1.xVec,data1.yVec,data1.elevationMaskedDetrended)
ylabel('Distance (m)')
hcb=colorbar('northoutside')
title(hcb, 'Detrended Elevation (m)')
colormap jet
caxis([-0.005 0.005])
set(gca, 'fontsize', 15)
nexttile
imagesc(data2.xVec,data2.yVec,data2.elevationMaskedDetrended)
ylabel('Distance (m)')
colormap jet
caxis([-0.005 0.005])
set(gca, 'fontsize', 15)
nexttile
imagesc(data3.xVec,data3.yVec,data3.elevationMaskedDetrended)
xlabel('Distance (m)')
ylabel('Distance (m)')
text(x,y,'t=9 hours', 'Color', 'w', 'Fontsize', 15)
colormap jet
caxis([-0.005 0.005])
set(gca, 'fontsize', 15)
nexttile
imagesc(data4.xVec,data4.yVec,data4.elevationMaskedDetrended)
xlabel('Distance (m)')
ylabel('Distance (m)')
colormap jet
caxis([-0.005 0.005])
set(gca, 'fontsize', 15)


