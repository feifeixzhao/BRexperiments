%% Load some of the data from Run 1
clear all; clc

[file_list, path_n]= uigetfile('.mat', 'grab the files', 'MultiSelect', 'on')

% get the file paths
if iscell(file_list) == 0 
    file_list = (file_list); 
end

for i = 1:length(file_list)
    filename=file_list(i);
    filename=char(filename);
    data_in(i)=load([path_n, filename], 'data');
    
    % loading x,y,and elevation data for selected times
    loc_x{i}=data_in(i).data.xVec; 
    loc_y{i}=data_in(i).data.yVec;
    loc_z{i}=data_in(i).data.elevationMasked;
    loc_zd{i}=data_in(i).data.elevationMaskedDetrended;
    loc_zraw{i}=data_in(i).data.elevationRaw;  
    timeHours{i}=data_in(i).data.metadata.timeHours;
end


%% drawing cross stream sections 
interval=11841; % 20m, X position of basin, the total length is 18528 (37 m)
[channelbelt1,channelbelt2,horz,horz_raw,lid,Zlid] = scourtime(interval,loc_z,loc_zraw);

% plotting evolution of channel
figure 
plot(movmean(horz_raw{1},60), 'Color', [0 0 0], 'LineWidth',2); 
hold on 
%yline(lid, '--', 'averaged floodplain elevation', 'LineWidth', 2, 'fontsize', 20);
plot(movmean(horz_raw{2},60), '--','Color', [0.4 0.4 0.4], 'LineWidth',2)
plot(movmean(horz_raw{3},60),'Color', [0.8 0.8 0.8], 'LineWidth',2)
%Zavg=nanmean(Zlid{1});
legend('t=1 hr', 't=6 hrs', 't=31 hrs')
set(gca, 'fontsize', 20)
xlabel('Cross Stream Profile')

% for i=1:length(horz) %plotting profiles 
%     figure(1)
%     hold on
%     plot(horz_raw{i},'Color', col(length(horz),:));
%     yline(lid,'-','channel datum');
%     set(gca, 'fontsize', 15)
%     xlabel('Cross Stream Position (m)')
%     ylabel('Raw Elevation (m)')
%     legendStrings = "t = " + string(timeHours)+ " hours";
%     legend(legendStrings, 'Location', 'Southeast')
% end

% plotting scour depths
for i=1:length(Zlid)
    figure(2)
    hold on
    %plotting probability density distributions 
    ksdensity(Zlid{i},'function','pdf')
    set(gca, 'fontsize',15)
    xlabel('Scour Depth')
    ylabel('Density')
end

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
%% analyzing space isolating time

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

% plotting all scour depths
% figure(1)
% ksdensity(Zlid(:))
% xlabel('Scour Depth')
% ylabel('Density')
% set(gca, 'fontsize', 15)
% 
% figure(2)
% imagesc(loc_x, loc_y, Zlid75)
% colorbar
% xlabel('Distance (m)')
% ylabel('Distance (m)')

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
