% NAME :
%           Select_Central_Coils(coils, coil_select_flag, verbose)
%
% DESCRIPTION:
%           Computes coil sensitivity
%
% INPUTS:
%           complex         coils                       coil sensitivity
%           logical         coil_select_flag            flag for selecting coils
%
% OUTPUTS:
%           complex         coils                       chosen coil sensitivity
%           double          sort_coil_ind               index for chosen coils
%
% PSEUDOCODE:
%           looks at central region of coil
%           choses top 50%
%
% NOTES:
%
%
% CHANGE LOG:
%     VERSION        DATE             AUTHOR                    CHANGE
%       1.1         2020-11-05       Datta Singh Goolaub        Creation



% Reconstruction code for radial PC MRI using Compressed Sensing, Motion
% Correction and Metric optimized Gating
%
% Datta Singh Goolaub (2020)
% University of Toronto / The Hospital For Sick Children

function [coils, sort_coil_ind] = Select_Central_Coils(coils, coil_select_flag, verbose,logf)
switch coil_select_flag
    case true
        % looks at max signal in the central portion of coil
        msk_coil_central_region = zeros(size(coils(:,:,1)));
        msk_coil_central_region(ceil(size(coils,1)/3):2*floor(size(coils,1)/3),ceil(size(coils,1)/3):2*floor(size(coils,1)/3)) = 1;
        msk_coil_central_region = repmat(msk_coil_central_region,[1 1 size(coils,3)]);
        msk_coil_central_region = msk_coil_central_region.*coils;
        msk_coil_central_region = squeeze(sum(sum(msk_coil_central_region,1),2));
        [~, sort_coil_ind] = sort(abs(msk_coil_central_region));
        sort_coil_ind = sort_coil_ind(round(0.5*size(coils,3)):end);
        coils = coils(:,:,sort_coil_ind);
        coils = (coils)./repmat(sqrt(sum(abs(coils).^2,3)),[1,1,size(coils,3)]);
        
        % displays chosen coils
        FETAL_LOGF (logf, verbose, 'Used coil:  %1.f.\n', sort_coil_ind)
    otherwise
        % index of all coils
        sort_coil_ind = 1:size(coils,3);
        
        % displays chosen coils
        FETAL_LOGF (logf, verbose, 'Using all coils.\n');

end