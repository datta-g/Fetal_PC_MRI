% NAME :
%           Sort_CINE_spokes(KSpace,traj,Times, vel_tag, RWaves, nFrames)
%
% DESCRIPTION:
%           Sorts data into CINE frames given a heart rate model
%
% INPUTS:
%           complex single  KSpace          matrix holding acquired data (samples x coils x radial views)
%           complex double  traj            trajectory corresponding to KSpace data
%           double          Times           vector with timestamps for each spoke (ms)
%           double          vel_tag         vector with velocity tags denoting ecode direction for each spoke
%           double          RWaves          vector with RR intervals (ms)
%           double          nFrames         number of frames in CINE
%
% OUTPUTS:
%           struct          rKpc            sorted k-space
%           struct          rk              sorted trajectory
%           struct          spokes_index    sorted indices
%
% PSEUDOCODE:
%           compute cardiac phases for each RR interval
%           assign data based on timestamp to each phase
%
% NOTES:
%
%
% CHANGE LOG:
%     VERSION        DATE             AUTHOR                    CHANGE
%       1.0         2020-11-05       Datta Singh Goolaub        Creation



% Reconstruction code for radial PC MRI using Compressed Sensing, Motion
% Correction and Metric Optimised Gating
%
% Datta Singh Goolaub (2020)
% University of Toronto / The Hospital For Sick Children

function [rKpc, rk, spokes_index] = Sort_CINE_spokes(KSpace,traj,Times, vel_tag, RWaves, nFrames)

% fills 1 missing beat
if (max(RWaves)<max(Times))
    RWaves(end+1) = RWaves(end)+(RWaves(end)-RWaves(end-1));
end

% prepares containers
for j = 1:size(traj,3)
    for i = 1:nFrames
        rKpc{i,j}=[];
        rk{i,j} = [];
    end
end

% loops over heart rate model and find cardiac phase for data
for i = 2:length(RWaves)
    minRtime = RWaves(i-1);
    maxRtime = RWaves(i);
    
    if ~((isnan(minRtime)) || (isnan(maxRtime)))
        edGes = linspace(minRtime,maxRtime,nFrames+1);
        Y(i-1,:) = discretize(Times,edGes);
    else
        Y(i-1,:)= repmat(NaN,[length(Times) 1]);
    end
end
Y(isnan(Y))=0;
Y = sum(Y,1);
Y(Y== nFrames+1) = nFrames;
[nY,oY] = sort(Y);
vel_tag=vel_tag(oY); %resorting v according to above

% assigns sorted data to output structs
for i = 1:nFrames
    rKpc{i,1}=KSpace(:,oY(nY==i & vel_tag == 1),:); % encode 1
    rk{i,1} = traj(:,oY(nY==i & vel_tag == 1)); % encode 1
    rk{i,2} = traj(:,oY(nY==i & vel_tag == 1)); % encode 2
    spokes_index{i,1} = oY(nY==i & vel_tag == 1);
    rKpc{i,2}=KSpace(:,spokes_index{i,1}+1,:); % encode 2
end
