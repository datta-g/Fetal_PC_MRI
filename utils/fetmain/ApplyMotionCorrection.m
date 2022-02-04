% NAME :  
%           ApplyMotionCorrection(KSpace, traj, Reg_transforms, nframes_RT, quiescent_period_range, PAR)
% 
% DESCRIPTION:
%           Runs the motion correction  for radial phase contrast MRI  
%
% INPUTS:
%           double          RT                          real-time series      
%           struct          PAR                         holds the reconstruction parameters
%
% OUTPUTS:
%           double          Reg_transforms              translation transforms
%           double          quiescent_period_range      indices for real time frames in quiescent period
%           double          frame_ref                   index of reference frame
% 
% PSEUDOCODE:
%           prepare real time series
%           select region of interest
%           find reference frame
%           run registration framework
%           find outliers
%           
% NOTES:
%           
%
% CHANGE LOG:
%     VERSION        DATE             AUTHOR                    CHANGE                                              
%       1.1         2020-11-05       Datta Singh Goolaub        Creation



% Reconstruction code for radial PC MRI using Compressed Sensing, Motion
% Correction and Metric Optimised Gating
%
% Datta Singh Goolaub (2020)
% University of Toronto / The Hospital For Sick Children

function [Corrected_KSpace, Quiescent_traj] = ApplyMotionCorrection(KSpace, traj, Reg_transforms, nframes_RT, quiescent_period_range, PAR)

% check if correction is to be performed
FETAL_LOGF (PAR.logf, PAR.verbose, '--- APPLY MOCO MODULE --- start.\n')

if ~PAR.MOCORT.perform
    Corrected_KSpace = KSpace;
    Quiescent_traj =traj;
    FETAL_LOGF (PAR.logf, PAR.verbose, 'Correction is off since registration is off.\n')
    return
end

% displays motion correction status
FETAL_LOGF (PAR.logf, PAR.verbose, 'Starting motion correction for %s.\n', PAR.Fname ); tic

% displays preparation of motion correction
FETAL_LOGF (PAR.logf, PAR.verbose, 'Interpolating motion parameters to acquisition TR.\n' )

% interpolate from temporal resolution of real time to repetition time
interped_transforms(1,:) = interp1(PAR.MOCORT.segment:PAR.MOCORT.segment*2:-(PAR.MOCORT.segment) + nframes_RT*PAR.MOCORT.segment*2,Reg_transforms(1,:),1:size(KSpace,2),'linear','extrap');
interped_transforms(2,:) = interp1(PAR.MOCORT.segment:PAR.MOCORT.segment*2:-(PAR.MOCORT.segment) + nframes_RT*PAR.MOCORT.segment*2,Reg_transforms(2,:),1:size(KSpace,2),'linear','extrap');


% applying translational shifts in k-space
Corrected_KSpace = Translate_Correct_KSpace(KSpace,interped_transforms, traj);


% take longest quiescience period and continue analysis
lengths = arrayfun(@(x) size(quiescent_period_range(x).indx,2), 1:numel(quiescent_period_range));
longest_quiescent_period = quiescent_period_range(lengths == max(lengths)).indx; % frame index of longest period

% spoke index of longest period
spokes_quies = (1 + (min(longest_quiescent_period)-1)*PAR.MOCORT.segment*2) : (max(longest_quiescent_period))*PAR.MOCORT.segment*2;

% keeping trajectory and corrected kspace for longest quiescent period
Corrected_KSpace = Corrected_KSpace(:,spokes_quies,:);
Quiescent_traj = traj(:,spokes_quies);


% save results
FileName=['Motion_Corrected_Data_' PAR.Fname ];  Pathname = PAR.Pathname;    SaveStyle =PAR.SaveInterMedRes;
SAVE_INTERMEDIATE_RESULTS(Pathname, FileName, SaveStyle, Corrected_KSpace, Quiescent_traj, interped_transforms);

% displays end of correction
FETAL_LOGF (PAR.logf, PAR.verbose, 'Raw data from %s has been motion corrected.\n', PAR.Fname )
FETAL_LOGF (PAR.logf, PAR.verbose, 'Correction time:  %1.f s.\n\n', toc)

FETAL_LOGF (PAR.logf, PAR.verbose, '--- APPLY MOCO MODULE --- end.\n \n')

end

