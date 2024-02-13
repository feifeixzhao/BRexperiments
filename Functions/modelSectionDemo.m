%demonstration of Doug's model
addpath utilities
clear variables

%Exner shape parameters
A = 4.3;
B = 4.3;
tba = 0.07; %Ambient shear
%diffusivity
Dg = 0.025;
%MPM transport coef
m = 1;
n = 1.5;
%Porosity
p = 0.4;


%Invariant
l_x = 20;
dx = 0.1; % cell size
x = 0:dx:l_x; % x
tc = tand(34); % theta critical for slope relaxation
nt = 20000; % number of time steps
eta = 0.001*rand(size(x)); %inital roughened bed
eta = eta - mean(eta);
dt = 0.001; % time step
E = 1; % avalanche diffusion

%vertical component of climb
dndt = 5e-6; %dz/dt

%indexing for periodic
x_index = 1:length(eta);
x_p1_index = [2:length(eta) 1];
x_m1_index = [length(eta) 1:length(eta)-1];

%storage array for topo-steps
Z = zeros(nt,length(eta));

%time loop for only 1D Douglas Jerolmack model
for idx = 1:nt
    
    %Exnerish shear stress - topography relation
    deta = (eta(x_index)-eta(x_m1_index))./dx;
    tx = tba.*(1+ A.*(eta-mean(eta)) + B.*deta);
    
    %shadow zone (lee faces)
    tx(tx < 0) = 0;

    %Avalanche
    AVcriteria = (eta(x_index)-eta(x_p1_index))./dx; 
    qa_x = E*(AVcriteria.^2 - tc^2).*AVcriteria.*(AVcriteria > tc);
    
    %MPM transport
    q =  m*(tx).^n + qa_x; %+ 0.2*rand(size(tx)); %randomize the flux a lil

    %upwind transport and topographic diffusion
    a_eta = (-dt/((1-p)*dx)).*(q(x_index)-q(x_m1_index));
    d_eta = (dt*Dg/(2*dx)).*(eta(x_p1_index)+eta(x_m1_index) - 2.*eta(x_index));
    
    %update topography
    eta = eta + a_eta + d_eta + dndt;
    
    %store timestep
    Z(idx,:) = eta;
    
    %uncomment to plot surface over time
    %{ 
    if mod(j,1000) == 0 
    
        plot(x,eta)
        title(num2str(j))
        drawnow

    end
    %}
end

%% plot section
smpInterval = 100;
figure('Position',[100,100,1800,900]);
[S,Z]=plot_Section(Z,x,smpInterval);

set(gca,'FontSize',14)
xlabel('Horizontal distance')
ylabel('Vertical distance')

%% animate section

%find max and min of topography through time
clim_max = max(Z(:));
clim_min = min(Z(:));

vidObj = VideoWriter('dougsModelAnimation','MPEG-4');
vidObj.Quality = 100;
% vidObj.Height = 1080;
% vidObj.Width = 1920;

open(vidObj)
figure('Position',[100,100,1800,900]);
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','Zbuffer');

idxStart = 100;
idxInterval = 100;

for idx=idxStart:idxInterval:nt-1
 

    [hlyr,hsrf]=section(topo2strat(Z(1:smpInterval:(idx+1),:)),[],x);

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
    %     pause(1e0)
end

%close the object and figure
close(vidObj);
close all
%% post-processing set-up
smpInt = 200; %sample intermittency
srtIdx = 300; %starting timestep to sample
idx2smp = srtIdx:smpInt:size(Z,1); %make it a vector

%% mean dune scales as a function of model time (modify function to calculate stdev, or any other value)
[mWaveLen,mHeight,mTrough,mCrest]=mScales(Z,dx,idx2smp);
eqHeight= eqScale(idx2smp.*dt,mHeight); %equilibrium height (could be any scale)

[mPeriod,WaveLen,Cel] = waveLenCel2(Z,idx2smp,dx,dt,'dtw');% I think there's an issue with te 

%% plot morphological stuff

%Dune height
figure('position',[500 500 800 800])
    plot(idx2smp.*dt,mHeight,'linewidth',2,'color','k')
    title('Dune height')
    legend({'dune height'})
    xlabel('time')
    ylabel('dune height')
    set(gca,'fontsize',12)
    
%Dune wavelength
figure('position',[500 500 800 800])
    plot(idx2smp.*dt,mWaveLen,'linewidth',2,'color','k')
    title('Dune wavelength')
    legend({'dune wavelength'})
    xlabel('time')
    ylabel('dune wavelength')
    set(gca,'fontsize',12)
    
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

%Set thickness through time
figure('position',[500 500 800 800])
    plot(idx2smp.*dt,cellfun(@expMean,R.setThickness),'linewidth',2,'color','k')
    title('Set thickness')
    legend({'{\mu_{setThickness}}'})
    xlabel('Time')
    ylabel('{\mu_{setThickness}}')
    set(gca,'fontsize',12)

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



