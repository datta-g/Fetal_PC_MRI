% NAME :
%           MOGRT (KSpace, Times, PAR)
%
% DESCRIPTION:
%           Real time reconstruction with compressed sensing
%
% INPUTS:
%           complex single  KSpace          matrix holding acquired data (samples x coils x radial views)
%           double          PAR             holds the reconstruction parameters
%           complex double  traj            trajectory corresponding to KSpace data
%
% OUTPUTS:
%           double          RT              real-time series for both flow encodes
%
% PSEUDOCODE:
%           compute coil sensitivity
%           perform initial reconstruction for flow comp
%           refine reconstruction with compressed sensing for flow comp
%           perform initial reconstruction for through plane encode
%           refine reconstruction with compressed sensing for through plane encode
%
% NOTES:
%
%
% CHANGE LOG:
%     VERSION        DATE             AUTHOR                    CHANGE
%       1.0         2020-11-05       Datta Singh Goolaub        Creation



% Reconstruction code for radial PC MRI using Compressed Sensing, Motion
% Correction and Metric Optimised Gating
%
% Datta Singh Goolaub (2020)
% University of Toronto / The Hospital For Sick Children

function RT = MOGRT(KSpace, traj, PAR)

% check if correction is to be performed
FETAL_LOGF (PAR.logf, PAR.verbose, '--- MOG RT RECON MODULE --- start.\n')
if ~PAR.MOGRT.perform
    RT = [];
    FETAL_LOGF (PAR.logf, PAR.verbose, 'MOG RT recon is off.\n')
    return
end

% storing number of readout samples
EchoLength = size(KSpace,1);
% creating a new kspace holder for coil sensitivity reconstruction
kspace_temp = KSpace(:,1:2:end,:);  traj_temp = traj(:,1:2:end);

% reshaping all data for reconstructions using gpuNUFFT
KSpace = reshape(KSpace,EchoLength*size(KSpace,2),size(KSpace,3),size(KSpace,4));
kspace_temp = reshape(kspace_temp,EchoLength*size(kspace_temp,2),size(kspace_temp,3),size(kspace_temp,4));

% displays start of reconstruction and coil sensitivity comptutation
FETAL_LOGF (PAR.logf, PAR.verbose, 'Starting high temporal resolution real time reconstruction for %s.\n', PAR.Fname)
FETAL_LOGF (PAR.logf, PAR.verbose, 'Computing coil sensitivity.\n'); tic;

% Computing coil sensitivity using Walsh et al approach
FT = gpuNUFFT_DSG([real(traj_temp(:)),imag(traj_temp(:))]',col(Weights_Radial_Data1(traj_temp,0)),2,3,8,[0.5*EchoLength,0.5*EchoLength],[],true);
STACKS = FT'*kspace_temp(:,:);STACKS = reshape(STACKS,[size(STACKS,1),size(STACKS,2),size(kspace_temp,2),size(kspace_temp,3)]);
coils = ismrm_estimate_csm_walsh(STACKS); % coil sensitivity

% Choosing coils with high central signal to lower load on computation
[coils, sort_coil_ind] = Select_Central_Coils(coils, PAR.MOGRT.coilsSelect, PAR.verbose, PAR.logf);
KSpace = KSpace(:,sort_coil_ind);
KSpace = KSpace./max(abs(KSpace(:)));

clear STACKS FT kspace_temp traj_temp %clearing unused variables

% Segmenting data into 2 real time series: one for each flow encode
KSpace = reshape(KSpace,[size(traj),size(KSpace,2)]);
[KSpace_resorted(:,:,:,:,1),traj_temp(:,:,:,1)] = Reconstruct_GRASP_RT_SW(KSpace(:,1:2:end,:),traj(:,1:2:end),(abs(traj(:,1:2:end))),PAR.MOGRT.segment,0);
[KSpace_resorted(:,:,:,:,2),traj_temp(:,:,:,2)] = Reconstruct_GRASP_RT_SW(KSpace(:,2:2:end,:),traj(:,2:2:end),(abs(traj(:,2:2:end))),PAR.MOGRT.segment,0);
KSpace = KSpace_resorted;  clear KSpace_resorted; %clearing unused variables
KSpace = permute(reshape(KSpace,[size(KSpace,1)*size(KSpace,2),size(KSpace,3),size(KSpace,4),size(KSpace,5)]),[1 2 5 3 4]);
w = Weights_Radial_Data1(traj_temp(:,:,:,1));

% displays gpu nufft set up
FETAL_LOGF (PAR.logf, PAR.verbose, 'Creating %d gpuNUFFT operators.\n', size(KSpace,4))
    
% Creating NUFFT operator for compressed sensing
for frame_in_realtime = 1:size(KSpace,4)
    param.GPU(frame_in_realtime) = gpuNUFFT_DSG([col(real(traj_temp(:,:,frame_in_realtime))),col(imag(traj_temp(:,:,frame_in_realtime)))]',col(w(:,:,frame_in_realtime)),2,3,8,[size(coils,1),size(coils,1)],(coils),true);
end

% start of recon
FETAL_LOGF (PAR.logf, PAR.verbose, 'Starting compressed sensing reconstruction for flow compensated data.\n')

% selecting data for real time series for flow compensated acquisition
param.y = KSpace(:,:,:,:,1);
% An initial reconstruction to initiate compressed sensing
RT_ini = param.GPU'*param.y;
% Transferring reconstruction parameters for compressed sensing
param.Weights = PAR.MOGRT.Weights*max(abs(RT_ini(:))); % scaling regularizers
param.Reg = PAR.MOGRT.REG;
param.nite = PAR.MOGRT.nite;
param.oneM = ones(size(param.y)); % mask
param.verbose = PAR.verbose;

% Compressed sensing reconstruction based on Lustig et al.
reset(parallel.gpu.GPUDevice.current()) % clear gpu to prevent overload
RT(:,:,:,1) = squeeze(CSL1NlCg(RT_ini,param));
% end of recon
FETAL_LOGF (PAR.logf, PAR.verbose, 'End of compressed sensing reconstruction for flow compensated data.\n')


% start of recon
FETAL_LOGF (PAR.logf, PAR.verbose, 'Starting compressed sensing reconstruction for through plane flow encode.\n')

% selecting data for real time series for flow compensated acquisition
param.y = KSpace(:,:,:,:,2);
% An initial reconstruction to initiate compressed sensing
RT_ini = param.GPU'*param.y;
param.Weights = PAR.MOGRT.Weights*max(abs(RT_ini(:))); % scaling regularizers
% Compressed sensing reconstruction based on Lustig et al.
reset(parallel.gpu.GPUDevice.current()) % clear gpu to prevent overload
RT(:,:,:,2) = squeeze(CSL1NlCg(RT_ini,param));
reset(parallel.gpu.GPUDevice.current()) % clear gpu to prevent overload
% end of recon
FETAL_LOGF (PAR.logf, PAR.verbose, 'End of compressed sensing reconstruction for through plane flow encode.\n')


% save results
FileName = ['RT_MOG_' PAR.Fname ];     Pathname = PAR.Pathname;    SaveStyle = PAR.SaveInterMedRes;
SAVE_INTERMEDIATE_RESULTS(Pathname, FileName, SaveStyle, RT);


% end of recon
FETAL_LOGF (PAR.logf, PAR.verbose, 'End of high temporal resolution real time reconstruction for %s has ended.\n', PAR.Fname)
FETAL_LOGF (PAR.logf, PAR.verbose, 'Reconstruction time:  %1.f s.\n\n', ceil(toc))
FETAL_LOGF (PAR.logf, PAR.verbose, '--- MOG RT RECON MODULE --- end.\n \n')

end
