% NAME :  
%           FETAL_HR_WITH_MOG(RT, Header, PAR)
% 
% DESCRIPTION:
%           Runs the metric optimized gating on real time reconstructions
%
% INPUTS:
%           double          RT                          real-time series      
%           struct          Header                      contains acqusition parameters
%           struct          PAR                         holds the reconstruction parameters
%
% OUTPUTS:
%           double          MOG_RWaveTimes              RR intervals (milliseconds)
% 
% PSEUDOCODE:
%           prepare real time series for metric optimized gating type
%           select region of interest using a dummy CINE
%           run multiparameter metric optimized gating algorithm
%           store heart rate model
%           
% NOTES:
%           to change the range of RR interval for seach, go to SET_RECON_PARAMS()           
%
% CHANGE LOG:
%     VERSION        DATE             AUTHOR                    CHANGE                                              
%       1.1         2020-11-05       Datta Singh Goolaub        Creation



% Reconstruction code for radial PC MRI using Compressed Sensing, Motion
% Correction and Metric optimized Gating
%
% Datta Singh Goolaub (2020)
% University of Toronto / The Hospital For Sick Children

function [MOG_RWaveTimes] = FETAL_HR_WITH_MOG(RT, Header, PAR)

% check if recon is to be performed
FETAL_LOGF (PAR.logf, PAR.verbose, '--- MOG MODULE --- start.\n')
if ~PAR.MOG.perform
    MOG_RWaveTimes = [];
    FETAL_LOGF (PAR.logf, PAR.verbose, 'MOG is off.\n')
    return
end

% displays MOG start
FETAL_LOGF (PAR.logf, PAR.verbose, 'Start of Metric Optimized Gating.\n'); tic;



% taking TR of each flow comp and encode acquisition
TR_2radialviews = Header{end}.hdr.MeasYaps.alTR{:}/1000;


% preparing data for ROI selection and MOG algorithm
dummy_fetal_hr = 385; % an initial heart rate only for ROI selection purposes
rt_duration = TR_2radialviews*PAR.MOGRT.segment; % temporal resolution of real time frames
time_series_frames = (rt_duration + TR_2radialviews)/2 : rt_duration : size(RT,3)*rt_duration - (rt_duration - TR_2radialviews)/2; % array of timestamps using TR of acqusition 
mog_input_images = PrepareDataForMOG(RT, PAR.MOG.METRIC); % prepare data according to set metric

% ROI selection query
FETAL_LOGF (PAR.logf, PAR.verbose, 'User selection for an ROI around target vessel.\n')

% a dummy CINE is created such that user can identify and select the vessel of interest
[ROI.y,ROI.x] = Select_ROI_CINE(ROI_Select_Idea((mog_input_images(:,:,1:end)),dummy_fetal_hr,max(time_series_frames),PAR.MOG.CardPhase,time_series_frames,1),{'CINE'});


% displays MOG start
FETAL_LOGF (PAR.logf, PAR.verbose, 'Searching for multiparameter heart rate model.\n')

% running metric optimized gating algorithm
[MOG_RWaveTimes,mog_resorted_CINE,log] = MRM_MOG_ISPACE(mog_input_images(ROI.y,ROI.x,:),time_series_frames,PAR,RT(45:130,45:120,:,:,:));

% displays details on heart rate
FETAL_LOGF (PAR.logf, PAR.verbose, 'Heart rate (mean +/- std) : %1.f +/- %1.f ms.\n', mean(diff(MOG_RWaveTimes)), std(diff(MOG_RWaveTimes)))



% save results
FileName=['MOG_FetalHeartRate_' PAR.Fname ];  Pathname = PAR.Pathname;    SaveStyle =PAR.SaveInterMedRes;
SAVE_INTERMEDIATE_RESULTS(Pathname, FileName, SaveStyle, MOG_RWaveTimes, mog_resorted_CINE, log, ROI);


% displays MOG end
FETAL_LOGF (PAR.logf, PAR.verbose, 'End of Metric Optimized Gating.\n')
FETAL_LOGF (PAR.logf, PAR.verbose, 'MOG time:  %1.f s.\n\n', toc)
FETAL_LOGF (PAR.logf, PAR.verbose, '--- MOG MODULE --- end.\n \n')

end