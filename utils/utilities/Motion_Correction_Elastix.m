% NAME :
%           Motion_Correction_Elastix(ImageRegistration,Mask,ElastixConfiguration,ElastixParam,ImageTarget,ImageMoving,frame_reference)
%
% DESCRIPTION:
%           Runs the registration on image frames
%
% INPUTS:
%           double          ImageRegistration           image series to be registered (x y time)
%           double          Mask                        mask for region used in registration
%           struct          ElastixConfiguration        configuration of elastix environment (see elxDefaultConfiguration)
%           struct          ElastixParam                elastix parameters
%           struct          ImageTarget                 reference image structure
%           struct          ImageMoving                 moving image structure
%           struct          frame_reference             reference frame index
%
% OUTPUTS:
%           int16           ImageRegistration           registered images
%           double          RegistrationTransforms      registration transforms
%           struct          RawTransforms               registration transforms output from elastix
%
% PSEUDOCODE:
%           run elastix module
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

function [ImageRegistration,RegistrationTransforms,RawTransforms] = Motion_Correction_Elastix(ImageRegistration,Mask,ElastixConfiguration,ElastixParam,ImageTarget,ImageMoving,frame_reference)

% prepare output transform struct
if strcmp(ElastixParam{1}.Transform,'EulerTransform')
    RegistrationTransforms=zeros(3,size(ImageRegistration,3));
elseif strcmp(ElastixParam{1}.Transform,'TranslationTransform')
    RegistrationTransforms=zeros(2,size(ImageRegistration,3));
end

% reference image
ImageTarget.Data = ImageRegistration(:,:,frame_reference);

Count=0;
for iFrame=1:size(ImageRegistration,3)
    if squeeze(sum(sum(ImageRegistration(:,:,iFrame))))~=0
        Count=Count+1;
        
        % moving image
        ImageMoving.Data=ImageRegistration(:,:,iFrame);
        % registration
        [IM, Transforms] = elxElastix(ElastixConfiguration,ElastixParam, ImageTarget, ImageMoving, 'FixedMask', Mask);
        RawTransforms{iFrame,1} = Transforms;
        
        % storing transforms
        if ~isempty(Transforms)
            RegistrationTransforms(:,iFrame)=Transforms{1,1}.TransformParameters;
            ImageRegistration(:,:,iFrame) = IM.Data;
        end
        
    end
end

RegistrationTransforms=RegistrationTransforms(:,1:Count);
end
