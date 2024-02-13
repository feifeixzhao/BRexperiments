% wrapper_prepare_dems_caesar.m
% Wrapper file for creating input DEMs for CAESAR.
% Created January 27, 2017 by Ajay Limaye (aslimaye@umn.edu).
% Last edited February 8, 2017 by Ajay Limaye (aslimaye@umn.edu).

addpath('source')

options.dx = 0.020677; % original DEM grid spacing
options.crop_lims.x = [0 36.42]; % crop original DEM using x1 x2 y1 y2, in meters. Do this, rather than row-column, to account for different DEM extents
options.crop_lims.y = [0.1 2.6];
options.downsample_factor = 10; % i.e., downsample by this factor
options.add_sidewalls = true;

infile='run4_processedElevation_010.mat';
outfile='caesar_10hr_10x';
prepare_dems_caesar(infile,outfile,options);


%% plot all distributions 

MaskedM_all=[MaskedM1;MaskedM2;MaskedM3;MaskedM4;MaskedM5];
figure(1)
subplot(2,2,1)
histogram(MaskedM_all,20)
title('Run 1')
xlabel('Flow depths (m)')
ylabel('Frequency')

subplot(2,2,2)
histogram(MaskedM_all,20)
title('Run 2')
xlabel('Flow depths (m)')
ylabel('Frequency')

subplot(2,2,3)
histogram(MaskedM_all,20)
title('Run 3')
xlabel('Flow depths (m)')
ylabel('Frequency')
ax = gca; % axes handle
ax.XAxis.Exponent = 0;


subplot(2,2,4)
histogram(MaskedM_all,20)
title('Run 4')
xlabel('Flow depths (m)')
ylabel('Frequency')

x =[MaskedM_1;MaskedM_2;MaskedM_3;MaskedM_4];
g = [zeros(length(MaskedM_1), 1); ones(length(MaskedM_2), 1); 2*ones(length(MaskedM_3), 1); 3*ones(length(MaskedM_4), 1)];
boxplot(x, g,'Labels',{'Run 1','Run 2','Run 3','Run 4'},'Symbol','')
ylabel('Flow Depth (m)')
ylim([0 0.015])

%% plot temporal change in distributions
x =[MaskedM_1;MaskedM_2;MaskedM_3;MaskedM_4;MaskedM_5];
g = [zeros(length(MaskedM_1), 1); ones(length(MaskedM_2), 1); 2*ones(length(MaskedM_3), 1); 3*ones(length(MaskedM_4), 1);4*ones(length(MaskedM_5), 1)];
figure(1)
subplot(2,2,1)
boxplot(x, g,'Labels',{'t = 10hr','t = 15 hr','t = 20 hr','t = 25 hr', 't = 30 hr'},'Symbol','')
ylabel('Flow Depth (m)')
title('Run 1')
ylim([0 0.015])

subplot(2,2,2)
boxplot(x, g,'Labels',{'t = 5hr','t = 6 hr','t = 7 hr','t = 8 hr', 't = 9 hr'},'Symbol','')
ylabel('Flow Depth (m)')
title('Run 2')
ylim([0 0.015])

subplot(2,2,3)
boxplot(x, g,'Labels',{'t = 5 hr','t = 10 hr','t = 15 hr','t = 20 hr', 't = 25 hr'},'Symbol','')
ylabel('Flow Depth (m)')
title('Run 3')
ylim([0 0.015])

subplot(2,2,4)
boxplot(x, g,'Labels',{'t = 3 hr','t = 7 hr','t = 11 hr','t = 15 hr', 't = 19 hr'},'Symbol','')
ylabel('Flow Depth (m)')
title('Run 4')
ylim([0 0.015])

%% plot flow depths temporally with dimensionless time 
% gather time steps
data0=load('Run1_WaterDepthMasked_1hr.mat')
data1=load('Run1_WaterDepthMasked_5hr.mat');
data2=load('Run1_WaterDepthMasked_7hr.mat');
data3=load('Run1_WaterDepthMasked_10hr.mat'); 
data4=load('Run1_WaterDepthMasked_15hr.mat'); 
data5=load('Run1_WaterDepthMasked_20hr.mat');

%concentate into matrix
Run1_medians=[median(data0.MaskedM);median(data1.MaskedM);median(data2.MaskedM);median(data3.MaskedM);median(data4.MaskedM);median(data5.MaskedM)];
Run4_quartile=[prctile(data1.MaskedM,95);prctile(data2.MaskedM,95);prctile(data3.MaskedM,95);prctile(data4.MaskedM,95);prctile(data5.MaskedM,95)];

%plotting
figure
plot(tstar, Run1_medians.*1000, '-o', 'color', [0.2 0.2 0.2],'MarkerSize', 10, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.2 0.2 0.2])
hold on
plot(tstar, Run2_medians.*1000, '-o', 'color', [0.8 0.8 0.8],'MarkerSize', 10, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.8 0.8 0.8])
plot(tstar, Run3_medians.*1000, '-o', 'color', [1 0 0],'MarkerSize', 10, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', [1 0 0])
plot(tstar, Run4_medians.*1000, '-o', 'color', [1 0.6 0.6],'MarkerSize', 10, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', [1 0.6 0.6])
plot(tstar,Run1_quartile*1000,'-s', 'color',[0.2 0.2 0.2], 'MarkerSize', 10, 'MarkerEdgeColor', 'black', 'MarkerFaceColor',[0.2 0.2 0.2])
plot(tstar,Run2_quartile*1000,'-s', 'color',[0.8 0.8 0.8], 'MarkerSize', 10, 'MarkerEdgeColor', 'black', 'MarkerFaceColor',[0.8 0.8 0.8])
plot(tstar,Run3_quartile*1000,'-s', 'color',[1 0 0], 'MarkerSize', 10, 'MarkerEdgeColor', 'black', 'MarkerFaceColor',[1 0 0])
plot(tstar,Run4_quartile*1000,'-s', 'color',[1 0.6 0.6], 'MarkerSize', 10, 'MarkerEdgeColor', 'black', 'MarkerFaceColor',[1 0.6 0.6])

%labels
xlabel('dimensionless time (t*)')
ylabel('flow depth (mm)')
xlim([1400 10600])
ylim([1.5 9])
legend({'Run 1 50th Percentile', 'Run 2 50th Percentile', 'Run 3 50th Percentile', 'Run 4 50th Percentile', 'Run 1 95th Percentile','Run 2 95th Percentile','Run 3 95th Percentile','Run 4 95th Percentile'},...
    'Location','eastoutside','NumColumns',2)
