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
%           double          resorted_CINE               resorted real-times into a CINE
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

function [MOG_RWaveTimes, resorted_CINE] = FETAL_HR_WITH_MOG(RT, Header, PAR)

% check if recon is to be performed
FETAL_LOGF (PAR.logf, PAR.verbose, '--- MOG MODULE --- start.\n')
if ~PAR.MOG.perform || PAR.PIPELINE.Pseudogating>0
    MOG_RWaveTimes = []; resorted_CINE = [];
    if PAR.PIPELINE.Pseudogating>0
        MOG_RWaveTimes = 0:PAR.PIPELINE.Pseudogating:(PAR.PIPELINE.Pseudogating+(Header{end}.hdr.Config.NLinMeas*Header{end}.hdr.MeasYaps.alTR{:}/1000));
    end
    FETAL_LOGF (PAR.logf, PAR.verbose, 'MOG is off.\n')
    return
end


for j = 1:size(RT,4)
    RT1(:,:,:,j) = imresize3(RT(:,:,:,j),[1/PAR.PIPELINE.ResolutionFraction*size(RT,1) 1/PAR.PIPELINE.ResolutionFraction*size(RT,2) size(RT,3)]);
end

% displays MOG start
FETAL_LOGF (PAR.logf, PAR.verbose, 'Start of Metric Optimized Gating.\n'); tic;



% taking TR of each flow comp and encode acquisition
TR_2radialviews = Header{end}.hdr.MeasYaps.alTR{:}/1000;


% preparing data for ROI selection and MOG algorithm
dummy_fetal_hr = 385; % an initial heart rate only for ROI selection purposes
rt_duration = TR_2radialviews*PAR.MOGRT.segment; % temporal resolution of real time frames
time_series_frames = (rt_duration + TR_2radialviews)/2 : rt_duration : size(RT1,3)*rt_duration - (rt_duration - TR_2radialviews)/2; % array of timestamps using TR of acqusition 
mog_input_images = PrepareDataForMOG(RT1, PAR.MOG.METRIC); % prepare data according to set metric

% ROI selection query
FETAL_LOGF (PAR.logf, PAR.verbose, 'User selection for an ROI around target vessel.\n')

% a dummy CINE is created such that user can identify and select the vessel of interest
[ROI.y,ROI.x] = Select_ROI_CINE(ROI_Select_Idea((mog_input_images(:,:,1:end)),dummy_fetal_hr,max(time_series_frames),PAR.MOG.CardPhase,time_series_frames,1),{'CINE'});


% displays MOG start
FETAL_LOGF (PAR.logf, PAR.verbose, 'Searching for multiparameter heart rate model.\n')

% running metric optimized gating algorithm
[MOG_RWaveTimes,mog_resorted_CINE,log] = MRM_MOG_ISPACE(mog_input_images(ROI.y,ROI.x,:),time_series_frames,PAR);

% displays details on heart rate
FETAL_LOGF (PAR.logf, PAR.verbose, 'Heart rate (mean +/- std) : %1.f +/- %1.f ms.\n', mean(diff(MOG_RWaveTimes)), std(diff(MOG_RWaveTimes)))

% resorts if needed
resorted_CINE = [];
if strcmp(PAR.PIPELINE.CINEType, 'resort')
    resorted_CINE(:,:,:,1) = resort_ISpaceRT(RT(:,:,:,1),time_series_frames,MOG_RWaveTimes,PAR.MOG.CardPhase);
    resorted_CINE(:,:,:,2) = resort_ISpaceRT(RT(:,:,:,2),time_series_frames,MOG_RWaveTimes,PAR.MOG.CardPhase);
end

% save results
FileName=['MOG_FetalHeartRate_' PAR.Fname ];  Pathname = PAR.Pathname;    SaveStyle =PAR.SaveInterMedRes;
SAVE_INTERMEDIATE_RESULTS(Pathname, FileName, SaveStyle, MOG_RWaveTimes, mog_resorted_CINE, log, ROI, resorted_CINE);


% displays MOG end
FETAL_LOGF (PAR.logf, PAR.verbose, 'End of Metric Optimized Gating.\n')
FETAL_LOGF (PAR.logf, PAR.verbose, 'MOG time:  %1.f s.\n\n', toc)
FETAL_LOGF (PAR.logf, PAR.verbose, '--- MOG MODULE --- end.\n \n')

end