% NAME :
%           Coil_Sens(rKpc,rk)
%
% DESCRIPTION:
%           Computes coil sensitivity
%
% INPUTS:
%           struct          rKpc                        kspace data
%           struct          rk                          trajectory
%
% OUTPUTS:
%           complex         coils                       coil sensitivity
%           complex         slices                      slice data from each coil
%
% PSEUDOCODE:
%           compute slice data
%           compute coil sensitivity with Walsh method
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

function [coils, slices] = Coil_Sens(rKpc,rk)

% Kspace and k prep
Ksp = [];
k = [];
for j = 1:size(rk,2)
    for i = 1:size(rk,1)
        k = cat(2,k,rk{i,j});
        Ksp = cat(2,Ksp,rKpc{i,j});
    end
end

% compute slice data from each coil
w = Weights_Radial_Data1(k,0);
EchoLength = size(Ksp,1);
Ksp=reshape(Ksp,EchoLength*size(Ksp,2),size(Ksp,3),size(Ksp,4));
FT = gpuNUFFT_DSG([real(k(:)),imag(k(:))]',w(:),2,3,8,[0.5*EchoLength,0.5*EchoLength],[],true);

% coil sensitivity using walsh
slices=FT'*Ksp(:,:);
slices=reshape(slices,[size(slices,1),size(slices,2),size(Ksp,2),size(Ksp,3)]);
[coils] = ismrm_estimate_csm_walsh(slices);
