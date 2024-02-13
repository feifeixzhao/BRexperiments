function [R] = fastPost(Zs,dt,idx2smp)
% fastPost2 - calculate bedform wavelength, height, celerity, peak and
% trough elevation through time, for chunks for time.

%results cell elements represent time chunk'ed populations for each simulation
%all time records
R.noDunes = {[]};
R.noSets = {[]};
%first timestep within currently preserved section
R.firstStep = {[]};
%all dune scales from firstStep to sampled timestep
R.duneHeights = {[]};
R.duneTroughs = {[]};
R.duneCrests = {[]};
%dune scales between firstSte
R.duneHeightsInt = {[]}; 
R.duneTroughsInt = {[]};
R.duneCrestsInt = {[]};
%stratigraphic scales 
R.setThickness = {[]};
R.stratTime = {[]};
R.surfTime = {[]};

[~,col] = size(Zs);

%iteration counter for storage
ic = 1;

%loop to calculate set thickness through time
 for kdx = 1:length(idx2smp)

    k = idx2smp(kdx);%end-timestep of sample
    timeVec = 1:k;
    
    %create kth timestep stratigraphy
    Z = topo2strat(Zs(1:k,:));
    Li = Z == Zs(1:k,:);
    Li_bs = logical_bounding_surface(Li,size(Z,1));
    [~, c] = max( Li ~=0, [], 1 ); %this finds the first timestep to become stratigraphy
    
    nstc = zeros(1,size(Z,2)); %number of sets
    stc = {[]}; %set thicknesses
    sfT = zeros(1,col);
    stT = zeros(1,col);
    
    
    %number of sets and set thickness for each node
    for i = 1:size(Z,2) %for each grid node
        sttp = Z(Li_bs(:,i),i); %find bounding surface elevatios
        bsElev = sttp(1:2:(end-2));
        stp = diff(bsElev); % this removes the set from the modern bedform
        stc{i} = stp(stp>0); % set thicknesses (non-uniform number)
        nstc(i) = numel(stc{i}); % number of sets per node
        
        %total time of strat and surfaces
        stT(i) = sum(Li(:,i)'.*timeVec>0).*dt; %preserved time
        sfT(i) = sum((~Li(:,i))'.*timeVec>0).*dt; %shreaded time
    end
    
    emtST = cellfun(@isempty,stc);
    stc(emtST) = [];
    
    R.setThickness{ic} = cat(1,stc{:}); % this records the set thicknesses for each investigation interval
    R.firstStep{ic}= c; % this records the first preserved strata for each column of Zs
    R.noSets{ic} = nstc;

    R.stratTime{ic} = stT;
    R.surfTime{ic} = sfT;
    ic = ic + 1;

 end

%mean periodicity, wavelength and celerity through time.
% [mPeriod,mWaveLen,Cel] = waveLenCel(Zs,smpInt,inputs.dx,inputs.dt,'dtw');

%find all peaks, troughs, bedform heights along domain for all chunks of
%model time

%save entire nodes geometry
htc = {[]}; %height
ctc = {[]}; %crests
ttc = {[]}; %trough

%save time of geometry
geoTime = {[]};


%loop through space, chunk in time later on!
 for j = 1:col
    %only consider topography from first bounding surface, onward. 
    Zsmp = Zs(:,j);
    [pks, plocs]=findpeaks(Zsmp, 'MinPeakDistance',10,'MinPeakProminence',.005);
    [tks, tlocs]=findpeaks(-Zsmp,'MinPeakDistance',10,'MinPeakProminence',.005);
    %logical index to make vectors same length
    li = 1:min(numel(tlocs),numel(plocs));

    %save entire nodes geometry
    htc{j} = pks(li) + tks(li); 
    ctc{j} = pks(li);
    ttc{j} = tks(li);
    
    %save time of geometry (choose to average times of peaks and troughs)
    geoTime{j} = mean([plocs(li),tlocs(li)],2);
 end
     

ic = 1;
    
%find values at set time intervals, aka, chunkify time for plotting
for kdx = 1:length(idx2smp)
    
    k = idx2smp(kdx); %kick out the index to sample
    
    %make logical index for dune geometry: first strata to k for each
    %chunk.
    fS = num2cell(R.firstStep{ic}); 
    
    %cell fun to chunkout, each cell should be first strata to end of time
    %interval.
    %unknown vector lengths, keep in cells.
    hts = cellfun(@(x,y,fs)(x((y>=fs) & (y<=k))),htc,geoTime,fS,'uni',0);
    cts = cellfun(@(x,y,fs)(x((y>=fs) & (y<=k))),ctc,geoTime,fS,'uni',0);
    tts = cellfun(@(x,y,fs)(x((y>=fs) & (y<=k))),ttc,geoTime,fS,'uni',0);
    %uniform output = 1
    ndtc = cellfun(@(x,y,fs)(numel(x((y>=fs) & (y<=k)))),htc,geoTime,fS);
    
    %sample crest heights over an interval of time.
    if kdx == 1
        htsInt = cellfun(@(x,y)(x(y<=k)),htc,geoTime,'uni',0);
        ctsInt = cellfun(@(x,y)(x(y<=k)),ctc,geoTime,'uni',0);
        ttsInt = cellfun(@(x,y)(x(y<=k)),ttc,geoTime,'uni',0);

    else
        htsInt = cellfun(@(x,y)(x((y>=idx2smp(kdx-1)) & (y<=k))),htc,geoTime,'uni',0);
        ctsInt = cellfun(@(x,y)(x((y>=idx2smp(kdx-1)) & (y<=k))),ctc,geoTime,'uni',0);
        ttsInt = cellfun(@(x,y)(x((y>=idx2smp(kdx-1)) & (y<=k))),ttc,geoTime,'uni',0);
      
    end
    
    
    
    %from last sample interval up-to current sample interval
    R.duneHeightsInt{ic} = cat(1,htsInt{:});
    R.duneCrestsInt{ic} = cat(1,ctsInt{:});
    R.duneTroughsInt{ic} = cat(1,ttsInt{:});

    %from first surface to sampled timestep
    R.duneHeights{ic} = cat(1,hts{:});
    R.duneCrests{ic} = cat(1,cts{:});
    R.duneTroughs{ic} = cat(1,tts{:});
    
    %number of dunes that has visited each grid node
    R.noDunes{ic} = ndtc;
        
    ic = ic + 1;
  
end