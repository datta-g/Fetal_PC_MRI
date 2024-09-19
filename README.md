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

<pre>
## Parameter file: SET_RECON_PARAMS.m
%% sets reconstruction parameters for reconstruction <br />
  PAR.SaveInterMedRes                 saves intermediate (0: not saving, 1: mat files) <br />
  PAR.verbose                         displays messages during reconstruction (0: off, 1: on) <br />
%% Data <br />
  PAR.dropSpokesTill                  drops initial spokes (must be even) <br />
%% Real time for motion correction <br />
  PAR.MOCORT.perform                  status 1 to perform RT recons, status 0 does not perform RT recon <br />
  PAR.MOCORT.regperform               status 1 to perform registration, status 0 does not perform registration <br />
  PAR.MOCORT.REG                      regularizer for compressed sensing <br />
  PAR.MOCORT.Weights                  coefficient for regularizer for compressed sensing <br />
  PAR.MOCORT.nite                     number of compressed sensing iterations <br />
  PAR.MOCORT.segment                  number of radial spokes in a real-time window <br />
  PAR.MOCORT.coilsSelect              flag to use use coils with high signal in central region <br />
  PAR.MOCORT.REGengine                'inbuilt' uses imreg functions from Matlab Library. 'elastix' uses elastix modules and needs installation <br />
%% Real time for metric optimized gating <br />
  PAR.MOGRT.perform                   status 1 to perform RT recons, status 0 does not perform RT recon <br />
  PAR.MOGRT.REG                       regularizer for compressed sensing <br />
  PAR.MOGRT.Weights                   coefficient for regularizer for compressed sensing <br />
  PAR.MOGRT.nite                      number of compressed sensing iterations <br />
  PAR.MOGRT.segment                   number of radial spokes in a real-time window <br />
  PAR.MOGRT.coilsSelect               flag to use use coils with high signal in central region <br />
%% metric optimized gating <br />
  PAR.MOG.perform                     status 1 to performs MOG, status 0 does not perform MOG <br />
  PAR.MOG.METRIC                      uses 'PC' setting for metric optimized gating <br />
  PAR.MOG.CardPhase                   number of cardiac phases in CINE <br />
  PAR.MOG.RRrange                     range of RR interval used to search fetal heart rate (milliseconds) <br />
%% CINE reconstruction <br />
  PAR.CINE.perform                    status 1 to perform CINE recons, status 0 does not perform CINE recon <br />
  PAR.CINE.REG                        regularizer for compressed sensing <br />
  PAR.CINE.nite                       number of compressed sensing iterations <br />
  PAR.CINE.CardPhase                  number of cardiac phases in CINE <br />
  PAR.CINE.coilsSelect                flag to use use coils with high signal in central region <br />
  PAR.CINE.Weights                    coefficient for regularizer 1 for compressed sensing <br />
%% Different reconstruction style <br />
  PAR.PIPELINE.AcquisFraction         Amount of data used as a fraction of overall acqusiition length [0.7 for an acquisition of 1000 spokes uses first 500 spokes in pipeline] <br />
  PAR.PIPELINE.ResolutionFraction     Reconstructed resolution ratio [0.5 is 50% of scanned resolution, max 1. min 0.1] <br />
  PAR.PIPELINE.MOGPARAM               'single' for single parameter MOG; 'multi' for multiparameter MOG <br />
  PAR.PIPELINE.CINEType               'resort' combines real-times into an estimate of CINE (quick not for analysis); 'CS' uses compressed sensing to compute CINE from raw data <br />
</pre>
