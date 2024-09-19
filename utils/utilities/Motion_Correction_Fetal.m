% NAME :
%           Motion_Correction_Fetal(Images,ROIy,ROIx,reference_frame)
%
% DESCRIPTION:
%           Runs the registration on image frames
%
% INPUTS:
%           double          Images                      image series to be registered (x y time)
%           double          ROIy                        vector with coordinates of edge of mask along y
%           double          ROIx                        vector with coordinates of edge of mask along x
%           double          reference_frame             reference frame index
%           str             regengine                   specification for registration modules ['inbuilt' for matlab library; 'elastix' for installed elastix software]
%
% OUTPUTS:
%           int16           ImageRegistration           registered images
%           double          RegistrationTransforms      registration transforms
%           struct          RawTransforms               registration transforms output from elastix % depreciated
%           double          ROIy                        vector with coordinates of edge of mask along y
%           double          ROIx                        vector with coordinates of edge of mask along x
%
% PSEUDOCODE:
%           prepare images series
%           imregister frames % depreciating elastix usage
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

function [ImageRegistration,RegistrationTransforms,RawTransforms,ROIy,ROIx] = Motion_Correction_Fetal(Images,ROIy,ROIx,reference_frame, regengine)

ImageRegistration=[];RegistrationTransforms=[];RawTransforms=[];

% ROI selection
if isempty(ROIy)
    [ROIy,ROIx] = Select_ROI_CINE(Images,{'CINE'});
end

% preparing images for elastix module
switch regengine
    case 'elastix'
        Images=127*abs(Images./abs(max(Images(:))));
        Images=int16(Images);
    case 'inbuilt'
        Images=255*abs(Images./abs(max(Images(:)))); % depreciating elastix usage
end

% imreg tools from matlab
switch regengine
    case 'elastix'
        % depreciating elastix usage for future pipelines
        % preparing configuration for registration
        ElastixConf = elxDefaultConfiguration;
        [ElastixParam,ImageTarget,ImageMoving,Mask] = mydefaultparameters(size(Images,1));
        Mask.Data=zeros(size(Images,1),size(Images,2),'int8');Mask.Data(ROIy,ROIx)=1;
        
        % cropping data to masked region
        for i=1:size(Images,3)
            itemp = Images(:,:,i);
            itemp(sum(Mask.Data,2)==0,:)=[];
            itemp(:,sum(Mask.Data,1)==0,:)=[];
            ImagesElastix(:,:,i) = itemp;
        end
        Mask.Data = ones(size(ImagesElastix,1),size(ImagesElastix,2));
        
        % running registration on data
        [ImageRegistration,RegistrationTransforms,RawTransforms] = Motion_Correction_Elastix(ImagesElastix,Mask,ElastixConf,ElastixParam,ImageTarget,ImageMoving,reference_frame);
        
    case 'inbuilt'
        RegistrationTransforms = CINE_Tool_MOCO( Images, repmat(Images(:,:,reference_frame),[1 1 size(Images, 3)]), ROIy,ROIx);
        ImageRegistration = CINE_Tool_Correct_RT(Images,-RegistrationTransforms,0);
        ImageRegistration = abs(ImageRegistration(ROIy,ROIx,:));
        RegistrationTransforms = -RegistrationTransforms';
end

end