% readme.m: Readme file for time series data for experimental braided channels, including
% example plotting scripts. Based on data from Limaye (2020), How do
% braided rivers grow channel belts?, JGR-Earth Surface,
% doi:10.1029/2020JF005570.
% Written by Ajay Limaye (ajay@virginia.edu), 12/23/20

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

% Load the 5th data file from run 1
load([pwd,'/run1/run1_processedElevation_002.mat'],'data')

% plot elevation masked by the extent of the channel belt
figure
imagesc(data.xVec,data.yVec,diff)
xlabel('Distance (m)')
ylabel('Distance (m)')
%axis image
h = colorbar;
ylabel(h, 'Elevation (m)')
colormap jet
caxis([-0.01 0.01])



%% to plot *detrended* elevation masked by the extent of the channel belt
figure
tiledlayout(3,1)
nexttile
imagesc(data1.xVec(9500:12018),data1.yVec(250:1000),data1.elevationMaskedDetrended(250:1000,9500:12018))
ylabel('Distance (m)')
x = max(data1.xVec(9500:12018))-1.3;% adapt the subtracted value to your needs
y = 0.6;
text(x,y,'t=3 hours', 'Color', 'w', 'Fontsize', 15)
hcb=colorbar('northoutside')
title(hcb, 'Detrended Elevation (m)')
colormap jet
caxis([-0.005 0.005])
set(gca, 'fontsize', 15)
nexttile
imagesc(data2.xVec(9500:12018),data2.yVec(250:1000),data2.elevationMaskedDetrended(250:1000,9500:12018))
ylabel('Distance (m)')
text(x,y,'t=6 hours', 'Color', 'w', 'Fontsize', 15)
colormap jet
caxis([-0.005 0.005])
set(gca, 'fontsize', 15)
nexttile
imagesc(data3.xVec(9500:12018),data3.yVec(250:1000),data3.elevationMaskedDetrended(250:1000,9500:12018))
xlabel('Distance (m)')
ylabel('Distance (m)')
text(x,y,'t=9 hours', 'Color', 'w', 'Fontsize', 15)
colormap jet
caxis([-0.005 0.005])
set(gca, 'fontsize', 15)
nexttile
imagesc(data4.xVec(9500:12018),data4.yVec(250:1000),data4.elevationMaskedDetrended(250:1000,9500:12018))
xlabel('Distance (m)')
ylabel('Distance (m)')
text(x,y,'t=25 hours', 'Color', 'w', 'Fontsize', 15)
colormap jet
caxis([-0.005 0.005])
set(gca, 'fontsize', 15)


