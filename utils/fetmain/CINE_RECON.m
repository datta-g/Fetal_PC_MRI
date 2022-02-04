% NAME :
%           CINE_RECON(Header, KSpace, traj, MOG_RWaveTimes, PAR)
%
% DESCRIPTION:
%           Reconstructs CINE for radial phase contrast MRI using compressed sensing
%
% INPUTS:
%           struct          Header          acquisition header
%           complex         KSpace          data matrix (read, ech, coil)
%           complex         traj            trajectory matrix
%           double          MOG_RWaveTimes  array of heart rate model
%           struct          PAR             reconstruction parameters
%
% OUTPUTS:
%           complex         CS              reconstructed CINE (x, y, time, velcity)
%
% PSEUDOCODE:
%           bins data into cardiac phases
%           compressed sensing for CINE reconstruction
%
%
% NOTES:
%           the values for the coefficients in this file were optimized using retrospectively undersampled data
%
% CHANGE LOG:
%     VERSION        DATE             AUTHOR                    CHANGE
%       1.0         2020-11-05       Datta Singh Goolaub        Creation



% Reconstruction code for radial PC MRI using Compressed Sensing, Motion
% Correction and Metric Optimised Gating
%
% Datta Singh Goolaub (2020)
% University of Toronto / The Hospital For Sick Children

function CS = CINE_RECON(Header, KSpace, traj, MOG_RWaveTimes, PAR)


% check if recon is to be performed
FETAL_LOGF (PAR.logf, PAR.verbose, '--- CINE RECON MODULE --- start.\n')
if ~PAR.CINE.perform
    CS = [];
    return
end

% displays motion correction status
FETAL_LOGF (PAR.logf, PAR.verbose, 'Starting CINE recon for %s.\n', PAR.Fname ); tic


% resort data
TR = Header{end}.hdr.MeasYaps.alTR{:}/2000; % grabing TR for an encode pair
acqusition_time = 0:TR:size(KSpace,2)*TR - TR; % building timestap array
v = repmat(1:2,[1 length(acqusition_time)/2]); % tags for encodes
[kspace_resorted_struct, traj_resorted_struct] = Sort_CINE_spokes(KSpace,traj,acqusition_time, v, MOG_RWaveTimes, PAR.CINE.CardPhase); clear KSpace v MOG_RWaveTimes;
[param.y,traj,~,max_r] = Change_struct_to_mat(kspace_resorted_struct,traj_resorted_struct);
% displays motion correction status
FETAL_LOGF (PAR.logf, PAR.verbose, 'Data resorted for CINE recon.\n')


% coil selection and sensitivity
[coils, ~] = Coil_Sens(kspace_resorted_struct,traj_resorted_struct); clear kspace_resorted_struct traj_resorted_struct;
[coils, sort_coil_ind] = Select_Central_Coils(coils, PAR.CINE.coilsSelect, PAR.verbose, PAR.logf);
param.y = param.y(:,1:max_r,sort_coil_ind,:,:);
param.y = param.y/max(abs(param.y(:)));
% displays motion correction status
FETAL_LOGF (PAR.logf, PAR.verbose, 'CINE Recon Coil sensitivity done.\n')
traj = traj(:,1:max_r,:,:,:);
%set up recon params
w = Weights_Radial_Data1(traj,0);
% w = abs(traj);
param.y=reshape(param.y,[size(param.y,1)*size(param.y,2),size(param.y,3),size(param.y,4),size(param.y,5)]);
param.oneM = w>0; % acquired data mask

param.oneM = reshape(param.oneM,size(param.oneM,1)*size(param.oneM,2),size(param.oneM,3),size(param.oneM,4));
param.oneM = repmat(param.oneM,[1 1 1 size(param.y,2)]);
param.oneM = permute(param.oneM,[1 4 5 2 3]);
% nufft operator set up
for time_dim=1:size(param.y,3)          %   time
    for vel_dim = 1:size(param.y,4)     %   velocity
        tempfft = gpuNUFFT_DSG([col(real(traj(:,:,time_dim))),col(imag(traj(:,:,time_dim)))]',col(w(:,:,time_dim)),2,3,8,[size(coils,1),size(coils,1)],(coils),true);
        param.GPU((time_dim-1)*size(w,4) + vel_dim) = tempfft;
    end
end
clear tempfft coils traj;
% displays motion correction status
FETAL_LOGF (PAR.logf, PAR.verbose, 'gpuNUFFT operators created.\n')


% reshape for nufft operations
param.y = permute(param.y,[1 2 5 3 4]);

% initial recon
CSre = param.GPU'*param.y;

% transferring CS parameters
param.Reg =[]; param.Weights=[];
param.Reg = PAR.CINE.REG;
param.Weights = PAR.CINE.Weights*max(abs(CSre(:)));
param.nite = PAR.CINE.nite;
param.verbose = PAR.verbose;

% displays motion correction status
FETAL_LOGF (PAR.logf, PAR.verbose, 'Starting compressed sensing.\n')

% iterative reconstruction
reset(parallel.gpu.GPUDevice.current())
CS = CSL1NlCg(CSre,param);
reset(parallel.gpu.GPUDevice.current())

% save results
FileName=['CINE_recon' PAR.Fname ];  Pathname = PAR.Pathname;    SaveStyle =PAR.SaveInterMedRes;
SAVE_INTERMEDIATE_RESULTS(Pathname, FileName, SaveStyle, CS);

% displays end of correction
FETAL_LOGF (PAR.logf, PAR.verbose, 'CINE series have been reconstructed.\n', PAR.Fname )
FETAL_LOGF (PAR.logf, PAR.verbose, 'Reconstruction time:  %1.f s.\n\n', toc)
FETAL_LOGF (PAR.logf, PAR.verbose, '--- CINE RECON MODULE --- end.\n \n')

end
