% NAME :
%           PrepareDataForMOG (data, mogtype)
%
% DESCRIPTION:
%           Prepares the real time images for metric optimized gating algorithm
%
% INPUTS:
%           double          data            real time series data (x y time velocity)
%           string          mogtype         type of preparation ('PC': phase contrast, 'CINE': anatomical)
%
% OUTPUTS:
%           double          MOG_prep_data   data prepared for metric optimized gating algorithm
%
% PSEUDOCODE:
%           -
%
% NOTES:
%           -
%
% CHANGE LOG:
%     VERSION        DATE             AUTHOR                    CHANGE
%       1.0         2020-11-05       Datta Singh Goolaub        Creation



% Reconstruction code for radial PC MRI using Compressed Sensing, Motion
% Correction and Metric Optimised Gating
%
% Datta Singh Goolaub (2020)
% University of Toronto / The Hospital For Sick Children

function MOG_prep_data = PrepareDataForMOG(data, mogtype)

switch mogtype
    % data is in the form of x-dim y-dim time vel
    case 'PC'
        % PC format
        data = (data(:,:,:,1).*conj(data(:,:,:,2)));
        MOG_prep_data = sqrt(abs(data)).*angle(data);
        
    case 'CINE'
        % CINE format
        data = (data(:,:,:,1).*conj(data(:,:,:,2)));
        MOG_prep_data = sqrt(abs(data));
end