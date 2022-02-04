% NAME :
%           SET_RECON_PARAMS ( )
%
% DESCRIPTION:
%           Sets reconstruction parameters for motion corrected radial phase contrast MRI
%
% INPUTS:
%           -
%
% OUTPUTS:
%           struct          PAR             holds the reconstruction parameters
%
% PSEUDOCODE:
%           assign parameter to field in PAR
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

function PAR = SET_RECON_PARAMS()

%% sets reconstruction parameters for reconstruction
PAR.SaveInterMedRes = 0;                % saves intermediate (0: not saving, 1: mat files)
PAR.verbose = 1;                        % displays messages during reconstruction (0: off, 1: on)

%% Data
PAR.dropSpokesTill = 0;                 % drops initial spokes (must be even)

%% Real time for motion correction
PAR.MOCORT.perform = 1;                 % status 1 to perform RT recons, status 0 does not perform RT recon
PAR.MOCORT.regperform = 1;              % status 1 to perform registration, status 0 does not perform registration
PAR.MOCORT.REG{1} = TVOP();             % regularizer 1 for compressed sensing
PAR.MOCORT.REG{2} = TV_Temp();          % regularizer 2 for compressed sensing
PAR.MOCORT.Weights(1) = 0.008;          % coefficient for regularizer 1 for compressed sensing
PAR.MOCORT.Weights(2) = 0.08;           % coefficient for regularizer 1 for compressed sensing
PAR.MOCORT.nite = 15;                   % number of compressed sensing iterations
PAR.MOCORT.segment = 64;                % number of radial spokes in a real-time window
PAR.MOCORT.coilsSelect = 1;             % flag to use use coils with high signal in central region

%% Real time for metric optimized gating
PAR.MOGRT.perform = 1;                  % status 1 to perform RT recons, status 0 does not perform RT recon
PAR.MOGRT.REG{1} = TVOP();              % regularizer 1 for compressed sensing
PAR.MOGRT.REG{2} = TV_Temp();           % regularizer 2 for compressed sensing
PAR.MOGRT.Weights(1) = 0.008;           % coefficient for regularizer 1 for compressed sensing
PAR.MOGRT.Weights(2) = 0.08;            % coefficient for regularizer 1 for compressed sensing
PAR.MOGRT.nite = 15;                    % number of compressed sensing iterations
PAR.MOGRT.segment = 8;                  % number of radial spokes in a real-time window
PAR.MOGRT.coilsSelect = 1;              % flag to use use coils with high signal in central region

%% metric optimized gating
PAR.MOG.perform = 1;                    % status 1 to performs MOG, status 0 does not perform MOG
PAR.MOG.METRIC = 'PC';                  % uses 'PC' setting for metric optimized gating
PAR.MOG.CardPhase = 15;                 % number of cardiac phases in CINE
PAR.MOG.RRrange = 350:550;              % range of RR interval used to search fetal heart rate (milliseconds)

%% CINE reconstruction
PAR.CINE.perform = 1;                   % status 1 to perform CINE recons, status 0 does not perform CINE recon
PAR.CINE.REG{1} = TVOP();               % regularizer 1 for compressed sensing
PAR.CINE.REG{2} = TV_Temp();            % regularizer 2 for compressed sensing
PAR.CINE.nite = 15;                     % number of compressed sensing iterations
PAR.CINE.CardPhase = PAR.MOG.CardPhase; % number of cardiac phases in CINE
PAR.CINE.coilsSelect = 1;               % flag to use use coils with high signal in central region
PAR.CINE.Weights(1) = 0.025;            % coefficient for regularizer 1 for compressed sensing
PAR.CINE.Weights(2) = 0.01;             % coefficient for regularizer 1 for compressed sensing

end