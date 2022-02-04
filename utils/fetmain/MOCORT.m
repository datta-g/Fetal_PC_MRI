% NAME :
%           MOCORT (KSpace, Times, PAR)
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
%           double          RT              real-time series
%
% PSEUDOCODE:
%           compute coil sensitivity
%           perform initial reconstruction
%           refine reconstruction with compressed sensing
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

function RT = MOCORT(KSpace, traj, PAR)

% check if recon is to be performed
FETAL_LOGF (PAR.logf, PAR.verbose, '--- MOCO RT RECON MODULE --- start.\n')
if ~PAR.MOCORT.perform
    RT = [];
    FETAL_LOGF (PAR.logf, PAR.verbose, 'MOCO RT recon is off.\n')
    return
end

KSpace = KSpace(:,1:2:end,:); % taking alternate data acquisitions (one encode)
traj = traj(:,1:2:end); % taking alternate data acquisitions (one encode)
EchoLength = size(KSpace,1); % number of samples in one readout

% displays start of reconstruction and coil sensitivity comptutation
FETAL_LOGF (PAR.logf, PAR.verbose, 'Starting low temporal resolution real time reconstruction for %s.\n', PAR.Fname)
FETAL_LOGF (PAR.logf, PAR.verbose, 'Computing coil sensitivity in MOCO RT recon.\n'); tic;

% Computing coil sensitivity using Walsh et al approach
KSpace = reshape(KSpace,EchoLength*size(KSpace,2),size(KSpace,3),size(KSpace,4));
FT = gpuNUFFT_DSG([real(traj(:)),imag(traj(:))]',col(Weights_Radial_Data1(traj,0)),2,3,8,[0.5*EchoLength,0.5*EchoLength],[],true);

STACKS = FT'*KSpace(:,:); STACKS = reshape(STACKS,[size(STACKS,1),size(STACKS,2),size(KSpace,2),size(KSpace,3)]);
coils = ismrm_estimate_csm_walsh(STACKS); % coil sensitivity

% Choosing coils with high central signal to lower load on computation
[coils, sort_coil_ind] = Select_Central_Coils(coils, PAR.MOCORT.coilsSelect, PAR.verbose, PAR.logf);
KSpace = KSpace(:,sort_coil_ind); KSpace=KSpace./max(abs(KSpace(:)));

clear STACKS FT %clearing unused variables

% Segmenting data into real time series
[KSpace,traj,~] = Reconstruct_GRASP_RT_SW(reshape(KSpace,[EchoLength,size(KSpace,1)/EchoLength,size(KSpace,2)]),traj,(abs(traj)),PAR.MOCORT.segment,0);
w = Weights_Radial_Data1(traj,0);

% displays gpu nufft set up
FETAL_LOGF (PAR.logf, PAR.verbose, 'Creating %d gpuNUFFT operators.\n', size(KSpace,4));

% Creating NUFFT operator for compressed sensing
for frame_in_realtime = 1:size(KSpace,4)
    param.GPU(frame_in_realtime) = gpuNUFFT_DSG([real(col(traj(:,:,frame_in_realtime))),imag(col(traj(:,:,frame_in_realtime)))]',col(w(:,:,frame_in_realtime)),2,3,8,[size(coils,1),size(coils,1)],coils,true);
end
KSpace = reshape(KSpace,[size(KSpace,1)*size(KSpace,2),size(KSpace,3),size(KSpace,4),size(KSpace,5)]);
KSpace = permute(KSpace,[1 2 5 3 4]);

% An initial reconstruction to initiate compressed sensing
RT = param.GPU'*KSpace;

% Transferring reconstruction parameters for compressed sensing
param.nite = PAR.MOCORT.nite;
param.Reg = PAR.MOCORT.REG;
param.Weights = PAR.MOCORT.Weights * max(abs(RT(:)));
param.y = KSpace;
param.oneM = ones(size(param.y)); % mask to indicate which field has data
param.verbose = PAR.verbose;

% displays start of recon
FETAL_LOGF (PAR.logf, PAR.verbose, 'Compressed sensing reconstruction has started. \n');

% Compressed sensing reconstruction based on Lustig et al.
reset(parallel.gpu.GPUDevice.current())     % clear gpu to prevent overload
RT = squeeze(CSL1NlCg(RT,param));           % compressed sensing reconstruction
reset(parallel.gpu.GPUDevice.current())     % clear gpu to prevent overload

% save results
FileName = ['RT_MOCO_' PAR.Fname ];     Pathname = PAR.Pathname;    SaveStyle = PAR.SaveInterMedRes;
SAVE_INTERMEDIATE_RESULTS(Pathname, FileName, SaveStyle, RT);

% end of recon
FETAL_LOGF (PAR.logf, PAR.verbose, 'End of low temporal resolution real time reconstruction for %s. \n', PAR.Fname);
FETAL_LOGF (PAR.logf, PAR.verbose, 'Reconstruction time:  %1.f s.\n\n', toc);
FETAL_LOGF (PAR.logf, PAR.verbose, '--- MOCO RT RECON MODULE --- end.\n \n')
