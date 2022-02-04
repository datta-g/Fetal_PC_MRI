% NAME :  
%           FETAL_RADIAL_PCMRI ( )
% 
% DESCRIPTION:
%           Runs the reconstruction for motion corrected radial phase contrast MRI    
% INPUTS:
%           - 
%
% OUTPUTS:
%           stored as mat files and dicom
% 
% PSEUDOCODE:
%           read raw data
%           prepare reconstruction parameters
%           reconstruct real-time series for motion correction
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

function FETAL_RADIAL_PCMRI()
% load dependencies
origPath = path;
addpath(genpath('./utils')); addpath(genpath('./class'))
resetPath = onCleanup(@()path(origPath));

% load parameters
PAR = SET_RECON_PARAMS();

% load data structure
[Header, KSpace, traj, ~, ~, PAR]=Read_Raw_Data([],PAR);
KSpace = KSpace(:,1:end-PAR.dropSpokesTill,:); traj = traj(:,1:end-PAR.dropSpokesTill);

%% Real time recon for motion correction: Compressed sensing with low temporal resolution data
RT = MOCORT(KSpace, traj, PAR);

%% Motion compensation: Extracting translational parameters
[Reg_transforms, quiescent_period_range] = RegisterRTFrames(RT, PAR);

%% Apply motion correction: Updating kspace and trajectory
[KSpace, traj] = ApplyMotionCorrection(KSpace, traj, Reg_transforms, size(RT,3), quiescent_period_range, PAR);

%% Real time recon for metric optimized gating: Compressed sensing with high temporal resolution data
RT = MOGRT(KSpace, traj, PAR);

%% Run metric optimized gating: Extracting heart rate from data
MOG_RWaveTimes = FETAL_HR_WITH_MOG(RT, Header, PAR);

%% Run CINE recon: Creating dynamic series for a representative heart beat
% MOG_RWaveTimes =0:500:20*500;
CS = CINE_RECON(Header, KSpace, traj, MOG_RWaveTimes, PAR);