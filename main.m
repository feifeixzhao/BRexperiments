% select DEMs to analyze
[file_list, path_n] = uigetfile('.mat', 'Select the files', 'MultiSelect', 'on');

% convert to cell array if not already
if ~iscell(file_list)
    file_list = {file_list}; 
end

% get the file paths
if iscell(file_list) == 0
file_list = (file_list);
end

%set up variables
for i = 1:length(file_list)
filename=file_list(i);
filename=char(filename);
data_in(i)=load([path_n, filename], 'data');
loc_x{i}=data_in(i).data.xVec; 
loc_y{i}=data_in(i).data.yVec;
loc_z{i}=data_in(i).data.elevationMasked;
loc_zd{i}=data_in(i).data.elevationMaskedDetrended;
loc_zraw{i}=data_in(i).data.elevationRaw;  
timeHours(i)=data_in(i).data.metadata.timeHours;
clean{i}=all(isnan(loc_z{i}),1);
loc_z{i}(:,clean{i})=[];
end

for j=1:length(file_list)
loc_zraw{j}(:,clean{j})=[];
loc_x{j}(clean{j})=[];
end

%% calculate scour depths
% Compute Zlid for each cell
Zlid = cellfun(@scourspace, loc_z, loc_zraw, loc_x, 'UniformOutput', false);

%% calculate set thickness
% setting up variables
time = timeHours(2:end);
xPosclip = loc_y{1}; 
x = abs(xPosclip-xPosclip(1));
dx = mean(diff(x));
dt = mean(diff(time));
nt=numel(time);

% calculate channel belt boundaries
for m=1:length(loc_z{1})
horz=cellfun(@(a) a(:,m), loc_z, 'UniformOutput', false);
horz_raw=cellfun(@(a) a(:,m), loc_zraw, 'UniformOutput', false);
channelbelt1=nan(1, size(loc_z, 2));
channelbelt2=channelbelt1;
    for i=1:length(horz) %plotting profiles 
        channelbelt1(i) = find(~isnan(horz{i}), 1, 'first');
        channelbelt2(i) = find(~isnan(horz{i}), 1, 'last');
    end
end

% processing data into cells
horz=cellfun(@cell2mat, horz, 'UniformOutput', false);
horz_raw=cellfun(@cell2mat, horz_raw, 'UniformOutput', false);
Z=cellfun(@(a) a.', horz_raw, 'UniformOutput', false);
Z=cellfun(@fliplr, Z, 'UniformOutput', false); 

% elevation values are indexed by last channel belt value
x1=cellfun(@(a) a(end,end), channelbelt1, 'UniformOutput', false); 
x2=cellfun(@(a) a(end,end), channelbelt2, 'UniformOutput', false); 
x1=cell2mat(x1);
x2=cell2mat(x2);

smpInterval = repmat({1},1,length(loc_z{1}));
x = repmat({x},1,length(loc_z{1}));
[S]=cellfun(@plot_Section, Z,x,smpInterval, 'UniformOutput', false); 

% calculating set thickness
S_new=cellfun(@(a) a(:,x1:x2), S, 'UniformOutput', false); 
diff_S=cellfun(@diff, diff_S, 'UniformOutput', false); 
loc_y=loc_y{1,1}(max(x1):max(x2)); % to plot the correct values for y axis
sthickness=cell2mat(diff_S); % unsorted set thicknesses
sthickness(sthickness <= 0.001) = []; % filtering out thicknesses from noise

%% plot synthetic strat section
interval=11841; % 20m, X position of basin, the total length is 18528 (37 m)

[channelbelt1, channelbelt2, horz, horz_raw] = scourtime(interval,loc_z,loc_zraw);

% converting to Roberts data demo 
%timeHours=cell2mat(timeHours');
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
% Z(1,:)=-1.445; % impose a lower boundary
Z_local=(Z-Z(1,1)).*1000; % convert to mm

% plot section
nt=numel(time);
smpInterval = 1;
figure(1)
subplot(4,1,4)
[S,Z_braid]=plot_Section(Z_local,x,smpInterval);
hold on
plot(x, Z_braid(5,:), 'k', 'LineWidth', 2)
xlim([0.98 1.71])
ylim([0 35])

mycolormap = customcolormap([0 .25 .5 .75 .97 1],{'#9d0142','#f66e45','#ffffbb','#65c0ae','#5e4f9f','#8f89a0'});
cb=colorbar();
colormap(mycolormap);

fig=figure(1)
h = axes(fig,'visible','off'); 
c = colorbar('north')
set(h,'Position',[0.35 0.92 0.3 0.065]);
c.YTick=c.Limits;
c.YTickLabel={'0', sprintf('10^%d', 4)}; %dimensionless time 
%% fit to Paola & Borgman model
setrnd = exprnd(3.0, 1, 1.9e6);
setrnd(setrnd<0.001)=[];

width = 0.001; % Vertical resolution
l = linspace(width, max(setrnd), 35);

% Calculate the PDF using Paola and Borgman's formula
a = 1.64493 / mean(setrnd);
f = @(s) ((a .* exp(-a * s)) .* (exp(-a * s) + (a * s) - 1)) ./ ((1 - exp(-a * s)).^2);
setpd = f(l);

% Create a bar graph for setrnd
[counts, edges] = histcounts(setrnd, 35); % Adjust the number of bins as needed

% Calculate the scaling factor for the PDF
scaling_factor = sum(counts) / trapz(l, setpd);

% Scale the PDF to match the bar graph
scaled_setpd = setpd * scaling_factor;

% Create the bar graph
figure
bar(edges(1:end-1), counts, 'BarWidth', 1, 'FaceColor', 'b');

% Plot the scaled PDF curve on the same plot
hold on;
plot(l, scaled_setpd, 'r', 'LineWidth', 2);
hold off;
xlim([-0.4 30])
xlabel('Value');
ylabel('Frequency (setrnd)');
title('Exponential Distribution Fit');

legend('setrnd', 'PDF (setpd)');

%% Fit an exponential distribution to sets
set=setthicknessall.*1000; % convert to mm
params = fitdist(set.', 'Exponential');

% Generate an x-axis for the exponential fit
x = linspace(1.5, max(set), 20);

% Calculate the corresponding exponential distribution values for the x-axis
y = pdf(params, x);

% Create a histogram (bar graph) of the data
figure;
histogram(set, 'NumBins', 20, 'Normalization', 'pdf', 'FaceAlpha', 0); % 'pdf' option to normalize the histogram


% Plot the exponential fit over the histogram
hold on; % This allows you to overlay the plot
plot(x, y, 'b', 'LineWidth', 1.7); % 'r' for red color, adjust line properties as needed

% Add labels and a legend
title('Run 4, Q = 0.12 L/s, S = 0.02');
legend('Observed sets', 'Predicted PDF');
ylabel('Probability density')
xlabel('set thickness,{\it s} (mm)')
xlim([0 30])


%% fitting set thicknesses to exponential

mean_set=3.1; % observed value
a_mean = 1.645/mean_set;
l = linspace(0, 25, 16);
f = @(s)((a_mean.*exp((-1).*a_mean.*s)).*(exp((-1).*a_mean.*s)+(a_mean.*s)-1))./((1-exp((-1).*a_mean.*s)).^2); % from Paola and Borgman
setPB=f(l); % pdf for predicted 

binCenters = l(1:end) - 0.5 * l(2);

set_obs=histcounts(setRun2.*1000, binCenters, 'Normalization', 'probability');
%set_obs= histcounts(setRun2.*1000, binCenters) / numel(setRun2);
figure;
% Plot the observed data as a bar graph with bin centers
bar(l, set_obs, 'FaceColor', 'red', 'DisplayName', 'Observed Data');
hold on
plot(l, setPB/sum(setPB), 'LineWidth', 2, 'Color', 'blue', 'DisplayName', 'Paola & Borgman PDF');

chi_2=sum((set_obs(2:end)-setPB(2:end)).^2./setPB(2:end)); % calculate chi-square
[h,p,stats] = chi2gof(l, 'Ctrs', l, 'Frequency', setpd_obs, 'Expected', setpd)

% calc goodness of fit
[h,p,stats] = chi2gof(l, 'Ctrs', l, 'Frequency', setpd_obs, 'Expected', setpd)
%[h,p,stats] = chi2gof(set1_braiding,'cdf',{@expcdf,mean(set1_braiding)}, 'NBins',30)
fprintf(['The chi squared value is ' num2str(stats.chi2stat) ' and the p value is ' num2str(p)])

% coefficient of variation
cv_sets=std(setRun4)/mean(setRun4);

%% fit scour depths to gamma
% calculating expected pdf 
Pd = fitdist(Zclean.','Gamma');
l=linspace(0, max(Zclean), 30);
gamma_pdf=pdf(Pd,l);

%calculating observed pdf
[heights, locations] = histcounts(Zlidall, 40);
width=0.001;
heights=heights/(length(Zlidall)*width); 

%plot
figure
% h=bar(locations, heights,'hist')
% set(h, 'FaceColor', 'w', 'EdgeColor', 'black');
plot(l, gamma_pdf)
hold on 
plot(locations,heights)
xlabel('Scour Depths (m)') 
ylabel('probability density')
set(gca, 'fontsize', 20)
legend('empirical PDF','gamma fit, alpha=2.004, beta=0.0029')
title('Run 4 Scour Depths')
saveas(gcf,'Run4Scour.fig')

% calc goodness of fit 
new_gam=gamrnd(Pd.a, Pd.b, 1000,1);
[h,p,stats] = chi2gof(new_gam,'cdf',{@gamcdf,Pd.a, Pd.b}, 'Nbins', 30);

figure
plot(l1*1000, P1./sum(l1),'LineWidth', 2, 'color', 'black')
hold on
plot(l2*1000, P2./sum(l2), 'LineWidth', 2, 'color', [0.5 0.5 0.5])
plot(l3*1000, P3./sum(l3), 'LineWidth',2, 'color', [1 0 0])
plot(l4*1000, P4./sum(l4), 'LineWidth',2,'color', [1 0.6 0.6])
legend('Run 1, Q = 0.25 L/s, S = 0.01', 'Run 2, Q = 0.12 L/s, S = 0.01', 'Run 3, Q = 0.25 L/s, S = 0.02', 'Run 4, Q = 0.12 L/s, S = 0.02')
xlabel('Scour depths (mm)')
ylabel('Probability density')
