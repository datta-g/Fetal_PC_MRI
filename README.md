# Fetal_PC_MRI
Reconstruction for radial PC MRI using Compressed Sensing, Motion Correction and Metric Optimized Gating


This set of code allows reconstructing accelerated fetal PC MRI data with motion correction and heart rate extraction steps.


Brief demo of usage is provided on https://sites.google.com/view/fetalpcmri/home

## Dependencies
GpuNufft      https://github.com/andyschwarzl/gpuNUFFT
(developed with V2.0.6)

Elastix       https://elastix.lumc.nl/
(developed with v4.800)


## Citation
Fetal Flow Quantification in Great Vessels Using Motion-Corrected Radial Phase Contrast MRI: Comparison With Cartesian. <br />
Goolaub et al.<br />
doi: 10.1002/jmri.27334


## Parameter file: SET_RECON_PARAMS.m
%% sets reconstruction parameters for reconstruction
  PAR.SaveInterMedRes                 saves intermediate (0: not saving, 1: mat files)
  PAR.verbose                         displays messages during reconstruction (0: off, 1: on)
%% Data
  PAR.dropSpokesTill                  drops initial spokes (must be even)
%% Real time for motion correction
  PAR.MOCORT.perform                  status 1 to perform RT recons, status 0 does not perform RT recon
  PAR.MOCORT.regperform               status 1 to perform registration, status 0 does not perform registration
  PAR.MOCORT.REG                      regularizer for compressed sensing
  PAR.MOCORT.Weights                  coefficient for regularizer for compressed sensing
  PAR.MOCORT.nite                     number of compressed sensing iterations
  PAR.MOCORT.segment                  number of radial spokes in a real-time window
  PAR.MOCORT.coilsSelect              flag to use use coils with high signal in central region
  PAR.MOCORT.REGengine                'inbuilt' uses imreg functions from Matlab Library. 'elastix' uses elastix modules and needs installation
%% Real time for metric optimized gating
  PAR.MOGRT.perform                   status 1 to perform RT recons, status 0 does not perform RT recon
  PAR.MOGRT.REG                       regularizer for compressed sensing
  PAR.MOGRT.Weights                   coefficient for regularizer for compressed sensing
  PAR.MOGRT.nite                      number of compressed sensing iterations
  PAR.MOGRT.segment                   number of radial spokes in a real-time window
  PAR.MOGRT.coilsSelect               flag to use use coils with high signal in central region
%% metric optimized gating
  PAR.MOG.perform                     status 1 to performs MOG, status 0 does not perform MOG
  PAR.MOG.METRIC                      uses 'PC' setting for metric optimized gating
  PAR.MOG.CardPhase                   number of cardiac phases in CINE
  PAR.MOG.RRrange                     range of RR interval used to search fetal heart rate (milliseconds)
%% CINE reconstruction
  PAR.CINE.perform                    status 1 to perform CINE recons, status 0 does not perform CINE recon
  PAR.CINE.REG                        regularizer for compressed sensing
  PAR.CINE.nite                       number of compressed sensing iterations
  PAR.CINE.CardPhase                  number of cardiac phases in CINE
  PAR.CINE.coilsSelect                flag to use use coils with high signal in central region
  PAR.CINE.Weights                    coefficient for regularizer 1 for compressed sensing
%% Different reconstruction style
  PAR.PIPELINE.AcquisFraction         Amount of data used as a fraction of overall acqusiition length [0.7 for an acquisition of 1000 spokes uses first 500 spokes in pipeline]
  PAR.PIPELINE.ResolutionFraction     Reconstructed resolution ratio [0.5 is 50% of scanned resolution, max 1. min 0.1]
  PAR.PIPELINE.MOGPARAM               'single' for single parameter MOG; 'multi' for multiparameter MOG
  PAR.PIPELINE.CINEType               'resort' combines real-times into an estimate of CINE (quick not for analysis); 'CS' uses compressed sensing to compute CINE from raw data
