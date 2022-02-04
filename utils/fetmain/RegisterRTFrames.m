% NAME :
%           RegisterRTFrames (RT, PAR)
%
% DESCRIPTION:
%           Runs the registration for radial phase contrast MRI
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

function [Reg_transforms, quiescent_period_range, frame_ref] = RegisterRTFrames(RT, PAR)

% check if reg is to be performed
FETAL_LOGF (PAR.logf, PAR.verbose, '--- REGISTRATION MODULE --- start.\n')
if ~PAR.MOCORT.regperform
    Reg_transforms = [];
    quiescent_period_range =[];
    frame_ref = [];
    FETAL_LOGF (PAR.logf, PAR.verbose, 'Registration is off.\n');
    return
end

% displays registration start
FETAL_LOGF (PAR.logf, PAR.verbose, 'Starting motion tracking pipeline for %s.\n', PAR.Fname )
tic

% preparing image and select region of interest
RT = abs(squeeze(RT));
for izz = 1:size(RT,3)
    RT(:,:,izz) = round(RT(:,:,izz)*255./max(max(RT(:,:,izz))));
end

% displays ROI selection query
FETAL_LOGF (PAR.logf, PAR.verbose, 'ROI selection for motion tracking for %s.\n', PAR.Fname )
    
% select region of interest
[ROI.y,ROI.x] = Select_ROI_CINE(Window_CINE(RT,0,0.75),{'CINE'});

% displays ROI selection end and registration start
FETAL_LOGF (PAR.logf, PAR.verbose, 'Successful ROI selection for motion tracking for %s.\n', PAR.Fname )
FETAL_LOGF (PAR.logf, PAR.verbose, 'Starting registration for %s.\n', PAR.Fname )


% frame of reference
frame_ref = find_reference_frame(RT, ROI);
% displays frame of reference
FETAL_LOGF (PAR.logf, PAR.verbose, 'Frame %d is frame of reference.\n', frame_ref )


% run registration
[RTR,Reg_transforms,allTransforms] = Motion_Correction_Fetal(RT,ROI.y,ROI.x,frame_ref);
[~,post_reg_cor,simi_mat] = find_reference_frame(double(RTR));
post_reg_cor = post_reg_cor/max(post_reg_cor);

% displays registration end
FETAL_LOGF (PAR.logf, PAR.verbose, 'Successful motion tracking for %s.\n', PAR.Fname )
FETAL_LOGF (PAR.logf, PAR.verbose, 'Locating outliers.\n')


% locate outliers
limitsCut = mean(post_reg_cor) - 1.5*iqr(post_reg_cor);
outlierIndex = post_reg_cor < limitsCut(1);


% displays registration end
FETAL_LOGF (PAR.logf, PAR.verbose, 'Frame ( %1.f ) were identified as outlier frames.\n',find(outlierIndex));
FETAL_LOGF (PAR.logf, PAR.verbose, 'Selecting quiescent period in scan.\n')


% dropping very bad motion here
% indices for frames in a long queiscent period are searched for
not_end_of_series = 1;
count_currnt_frame = 0;
quiescent_period_tag = 1;
quiescent_period_range(1).indx  = []; % creating a holder for indices
prev_status_tag = 1;
cut_off_for_accepted_period_length = 7; % 12 RT is equivalent to 760 spokes
while (not_end_of_series)
    count_currnt_frame = count_currnt_frame + 1;
    
    if outlierIndex(count_currnt_frame) == 0 && prev_status_tag == 1
        quiescent_period_range(quiescent_period_tag).indx = [quiescent_period_range(quiescent_period_tag).indx count_currnt_frame];
    elseif outlierIndex(count_currnt_frame) == 0
        quiescent_period_range(quiescent_period_tag).indx = [count_currnt_frame];
        prev_status_tag = 1;
    elseif prev_status_tag == 1
        if length(quiescent_period_range(quiescent_period_tag).indx) < cut_off_for_accepted_period_length
            quiescent_period_range(quiescent_period_tag).indx=[];
        elseif length(quiescent_period_range(quiescent_period_tag).indx) >= cut_off_for_accepted_period_length
            quiescent_period_tag = quiescent_period_tag + 1;
            quiescent_period_range(quiescent_period_tag).indx = [];
        else
            if count_currnt_frame < length(post_reg_cor)-cut_off_for_accepted_period_length+1
                quiescent_period_tag = quiescent_period_tag + 1;
                quiescent_period_range(quiescent_period_tag).indx = [];
            end
        end
        prev_status_tag = 0;
    end
    
    if count_currnt_frame == length(post_reg_cor)
        if length(quiescent_period_range(quiescent_period_tag).indx) < cut_off_for_accepted_period_length
            quiescent_period_range(quiescent_period_tag) = [];
        end
        not_end_of_series = 0;
    end
end


% save results
FileName=['MOCO_Details' PAR.Fname ];  Pathname = PAR.Pathname;    SaveStyle = PAR.SaveInterMedRes;
SAVE_INTERMEDIATE_RESULTS(Pathname, FileName, SaveStyle, RTR, allTransforms, Reg_transforms, frame_ref, ROI, quiescent_period_range, post_reg_cor, simi_mat, outlierIndex);


% displays registration end
FETAL_LOGF (PAR.logf, PAR.verbose, 'End of motion tracking pipeline for %s.\n', PAR.Fname )
FETAL_LOGF (PAR.logf, PAR.verbose, 'Motion tracking time:  %1.f s.\n\n', toc)
FETAL_LOGF (PAR.logf, PAR.verbose, '--- REGISTRATION MODULE --- end.\n \n')

end
