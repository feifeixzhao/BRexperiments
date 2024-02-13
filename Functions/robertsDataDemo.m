%Function to load Robert's data and visualize 
%load data
clear all
clc
close all

load('Limaye2020.mat');

time = timeHours(2:end);
x = abs(xPosclip-xPosclip(1));
Z = Depth;
Z = fliplr(Z); %comment out for working with Guala data
dx = mean(diff(x));
dt = mean(diff(time));
nt=numel(time);

%% plot section
nt=numel(time);
smpInterval = 1;
figure('Position',[100,100,1800,900])
S=plot_Section(Z,x,smpInterval);
set(gca,'FontSize',20)
xlabel('Horizontal distance')
ylabel('Vertical distance')

%%
figure
% subplot(2,1,1)
[F,hsrf,hlyr]=section(S(1:smpInterval:end,:),[],x);
subplot(2,1,2)
[F,hsrf]=sectioncol(S(1:smpInterval:end,:),[],x);

hold on

FX=[];
for i=1:size(S,2);
    FX(:,i)=central_diff(S(:,i),(0:dt:dt*(size(S,1)-1)));
end
%FX=FX';
Ss=S;
Ss(0<=FX)=nan;
Smsk=[];
Smsk=repmat(S(end,:),size(S,1),1);
Ss(Ss==Smsk)=nan;

plot(x,Ss,'color','k')
axis 'tight'

%% animate section
%find max and min of topography through time
clim_max = max(Z(:));
clim_min = min(Z(:));

vidObj = VideoWriter('robertsDataAnimation','MPEG-4');
vidObj.Quality = 100;
% vidObj.Height = 1080;
% vidObj.Width = 1920;


open(vidObj)
figure('Position',[100,100,1800,900]);
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','Zbuffer');

idxStart = 1;
idxInterval = 1;

for idx=idxStart:idxInterval:nt-1
 

    section(topo2strat(Z(1:smpInterval:(idx+1),:)),[],x);

    %fix aspect ratio if desired
    %daspect([1,1/300,1])
    ylim([min(Z(:)) max(Z(:))])
    title(['\Delta t  x ' num2str(idx)])
    set(gca,'FontSize',14)
 
    drawnow;
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
    cla
    
    %pause a little, if you want to watch as the video is made
   pause(1e0)
end

%close the object and figure
close(vidObj);
close all

%% post-processing set-up
smpInt = 1; %sample intermittency
srtIdx = 1; %starting timestep to sample
idx2smp = srtIdx:smpInt:size(Z,1); %make it a vector


%% mean dune scales as a function of model time (modify function to calculate stdev, or any other value)
[mWaveLen,mHeight,mTrough,mCrest]=mScales(Z,dx,idx2smp);
eqHeight= eqScale(idx2smp.*dt,mHeight); %equilibrium height (could be any scale)

%% plot morphological stuff

%Dune height
figure('position',[500 500 800 800])
    plot(idx2smp.*dt,mHeight,'linewidth',2,'color','k')
    title('Dune height')
    legend({'dune height'})
    xlabel('time')
    ylabel('dune height')
    set(gca,'fontsize',12)
    savefig('duneheight')
%Dune wavelength
figure('position',[500 500 800 800])
    plot(idx2smp.*dt,mWaveLen,'linewidth',2,'color','k')
    title('Dune wavelength')
    legend({'dune wavelength'})
    xlabel('time')
    ylabel('dune wavelength')
    set(gca,'fontsize',12)
    savefig('wavelenth')
%Dune topographic envelope
figure('position',[500 500 800 800])
    plot(idx2smp.*dt,mTrough,'linewidth',2,'color','k')
    hold on
    plot(idx2smp.*dt,mCrest,'linewidth',2,'color','r')
    title('Dune height')
    legend({'dune trough','dune crest'})
    xlabel('time')
    ylabel('dune crest')
    set(gca,'fontsize',12)
    savefig('topoenvelope')


%% calculate stratigraphic information (may take a minute)

%Guide: R - a 'R'esults structure has the following fields:
%
% Total results for stratigraphic section upto timestep of interest (each element in idx2smp)
%
%noDunes: the number of dunes that have visited each grid node, for each sample interval 
%noSets: the number of sets stacked vertically above each grid node, for each sample interval 
%firstStep: this is the timestep that corresponds to the lowest(earliest) bounding surface in each vertical section about each grid node
%duneHeights: all heights for time after 'firstStep' and up to the sample interval
%duneTroughs: all dune trough elevations for time after 'firstStep' and up to the sample interval
%duneCrests: all dune crest elevations for time after 'firstStep' and up to the sample interval
%setThickness: set thickness for time after 'firstStep' and up to the sample interval 
%stratTime: the preserved time in each vertical section
%surfTime: the time in SURFaces, or shredded time in each vertical section
%
% Chunk'ed results: gives a sense of change during a simulation/experiment
%
%duneHeightsInt: dune heights calculated for time between sample intervals
%duneTroughsInt: dune trough elevation calculated for time between sample intervals
%duneCrestsInt: dune crest elevation calculated for time between sample intervals

R = fastPost(Z,dt,idx2smp);
%% plot stratigraphic information (just examples)

%Fractional dune preservation as a function of time
figure('position',[500 500 800 800])
    plot(idx2smp.*dt,cellfun(@mean,R.noDunes),'linewidth',2,'color','k')
    hold on
    plot(idx2smp.*dt,cellfun(@mean,R.noSets),'linewidth',2,'color','r')
    title('Fractional dune preservation')
    legend({'average #dunes per node','average #sets per node'})
    xlabel('Time')
    ylabel('count')
    set(gca,'fontsize',12)
    savefig('preservation')
%Set thickness through time
figure('position',[500 500 800 800])
    plot(idx2smp.*dt,cellfun(@expMean,R.setThickness),'linewidth',2,'color','k')
    title('Set thickness')
    legend({'{\mu_{setThickness}}'})
    xlabel('Time')
    ylabel('{\mu_{setThickness}}')
    set(gca,'fontsize',12)
    savefig('setthickness')
%shredded and preserved time plotted against time
figure('position',[500 500 800 800])
    plot(idx2smp.*dt,cellfun(@mean,R.stratTime),'linewidth',2,'color','k')
    hold on
    plot(idx2smp.*dt,cellfun(@mean,R.surfTime),'linewidth',2,'color','r')
    title('Time preservation')
    legend({'Preserved time','Shredded time'})
    xlabel('Total Time')
    ylabel('Time Fraction')
    set(gca,'fontsize',12)
    savefig('shredding')
%Time ratio (preserved/shredded) as a function of time
figure('position',[500 500 800 800])
    v2p = cellfun(@mean,cellfun(@(x,y)(x./y),R.stratTime,R.surfTime,'uni',0));
    x2p = idx2smp.*dt;
    loglog(x2p,v2p,'linewidth',2,'color','k')
    hold on
    ft = fittype('power1');
    li = isfinite(v2p);
    fo = fit(x2p(li)',v2p(li)',ft);
    plot(x2p(li),fo(x2p(li)),'b--')
    title('Time preservation ratio')
    legend({'Time ratio', 'powerlaw Fit'})
    xlabel('Total Time')
    ylabel('Time ratio')
    set(gca,'fontsize',12)
    savefig('preservationratio')
%% Combine topographic and stratigraphic information for comparison to Paola & Borgman Theory or.. whatever

%for example, let's investigate variation in the reference level for Paola
%and Borgman (1991) and see if we can find a better match for Robert's data
%this is very similar to Brige and Leclair (2002)? I think?

%for reference, Paola and Borgman set the reference level to 2, indicating
%that dune topography below mean elevation contribute to set generation. 


figure('position',[500 500 800 800])
    omega = cellfun(@expMean,R.setThickness)./cellfun(@gammaMean,R.duneHeights); % mean set thickness from MLE exponential distribution / mean dune height from MLE gamma distribution
    cv = cellfun(@gammaCov,R.duneHeights); %coefficent of variation of dune height (MLE
    cv2 = cv.^2;
    refLevel = 2:6;
    loglog(cv2,omega,'ko')
    hold on
    xtp = linspace(min(cv2),max(cv2));
    for idx = 1:length(refLevel)
        loglog(xtp,(1.645/refLevel(idx)).*xtp,'k--')
        text(xtp(end-20),(1.645/refLevel(idx)).*xtp(end-20),['Ref Level = ',num2str(refLevel(idx))])
    end
 
   
    title('Stochastic Preservation Investigation')
    legend({'Data','Paola and Borgman (1991)'})
    xlabel('cv^{2}')
    ylabel('omega')
    set(gca,'fontsize',12)
    savefig('PaolaBorgman')