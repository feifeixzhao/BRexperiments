function [mWaveLen,mHeight,mTrough,mCrest]=mScales(Z,dx,idx2smp)
%function to calculate average values of bedform scales in space, for a
%moment in time. A "better" way to do this is shwon in fastPost.m, where each
%peak is found in a time-series at each node location.

    mWaveLen = zeros(size(idx2smp)); %preallocate
    mHeight = zeros(size(idx2smp));
    mTrough = zeros(size(idx2smp));
    mCrest = zeros(size(idx2smp));
    
    ic = 1;
    
    for jdx = 1:length(idx2smp) %loop over indices to sample

        idx = idx2smp(jdx); %kick out the index
        z = Z(idx,:); %truncate topography array

        [pks,plocs] = findpeaks(z); %find all peaks in space, for an instant in time
        mWaveLen(ic) = nanmean(diff(plocs))*dx; %calculate wavelength from peaks

        [tks,~] = findpeaks(-z);

        tks = -tks;

        mTrough(ic) = nanmean(tks); %find ar
        mCrest(ic) = nanmean(pks);
        mHeight(ic) = mCrest(ic) + abs(mTrough(ic));
        
        ic = ic + 1;
    end