% NAME :  
%           Change_struct_to_mat(rKpc,rk)
% 
% DESCRIPTION:
%           Transforms struc to mat
%
% INPUTS:
%           struct          rKpc            sorted k-space
%           struct          rk              sorted trajectory
%
% OUTPUTS:
%           complex         KSpace          matrix with kspace data (samples x profiles x coils x time x vel)
%           complex         traj            matrix with trajectory (samples x profiles x time x vel)
% 
% PSEUDOCODE:
%           
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


function [KSpace,traj,ws,min_r] = Change_struct_to_mat(rKpc,rk)

maxProf = 0;
for i = 1:size(rKpc,1)
    for j = 1:size(rKpc,2)
        if maxProf < size(rKpc{i,j},2)
            maxProf = size(rKpc{i,j},2);
        end
    end
end
KSpace = zeros(size(rKpc{1,1},1),maxProf,size(rKpc{1,1},3));
traj =  zeros(size(rKpc{1,1},1),maxProf);

min_r = Inf;
for i = 1:size(rKpc,1)
    for j = 1:size(rKpc,2)
        str = 1;
        ws(1,j,i) = str;
        KSpace(:,str:(str-1+size(rk{i,j},2)),:,i,j) = rKpc{i,j};
        traj(:,str:(str-1+size(rk{i,j},2)),i,j) = rk{i,j};
        str = str + size(rk{i,j},2);
        ws(2,j,i) = str-1;
        min_r = min(min_r,str-1);
    end
end