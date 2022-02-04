% NAME :
%           Read_Raw_Data (fileToRecon, PAR)
%
% DESCRIPTION:
%           Read radial phase contrast MRI raw data using mapVBVD
%
% INPUTS:
%           string          fileToRecon     contains filepath and filename. If empty (''), user browses to file
%           struct          PAR             contains parameters for recon
% OUTPUTS:
%           struct          Header          holds acquisition parameters and data description
%           double          Times           matrix holding timestamps for each acquisition
%           complex single  KSpace          matrix holding acquired data (samples x coils x radial views)
%           double          traj            trajectory matrix in the form kx + i*ky (samples x radial views)
%
% PSEUDOCODE:
%           -
%
% NOTES:
%           -
%
% CHANGE LOG:
%     VERSION        DATE             AUTHOR                    CHANGE
%       1.0         2020-11-05       Datta Singh Goolaub        Creation



% Reconstruction code for radial PC MRI using Compressed Sensing, Motion
% Correction and Metric Optimised Gating
%
% Datta Singh Goolaub (2020)
% University of Toronto / The Hospital For Sick Children

function [Header, KSpace, traj, Times, Noise, PAR] = Read_Raw_Data(fileToRecon, PAR)


if ~exist('fileToRecon','var') || isempty(fileToRecon)
    info = 'Select raw data file for reconstruction';
    [fname,pathname]=uigetfile('*.dat',info);
    fileToRecon = [pathname fname];
end


% prepare path and filename for saving outputs
[pathname, fname] = fileparts(fileToRecon);
PAR.Pathname = pathname;
PAR.Fname = fname;
PAR.logf = [fname '_log.txt'];
FETAL_LOGF (PAR.logf, PAR.verbose, '--- READ RAW DATA MODULE --- start.\n')

% start read
Header = mapVBVD(fileToRecon);

if length(Header) == 1
    H{1,1} = Header;
    Header = H;
end

for iMeas=1:length(Header)
    
    if isfield(Header{1,iMeas},'image')
        KSpace = Header{1,iMeas}.image.unsorted();
        KSpace=squeeze(KSpace);
        KSpace = permute(KSpace,[1 3 2]); % transform from   readout x coils x radial views   to    readout x radial views x coils
        for i=1:length(Header{1,iMeas}.image.Lin)
            Times(Header{1,iMeas}.image.Lin(i),Header{1,iMeas}.image.Set(i),Header{1,iMeas}.image.Sli(i),Header{1,iMeas}.image.Ave(i),Header{1,iMeas}.image.Rep(i))=2.5*Header{1,iMeas}.image.timestamp(i);
        end
        Times=squeeze(Times);
    else
        KSpace=[];
        Times=[];
    end
    
    if isfield(Header{1,iMeas},'noise')
        Noise = Header{1,iMeas}.noise();
        
        for i=1:length(Header{1,iMeas}.noise.Lin)
            Times(Header{1,iMeas}.noise.Lin(i),Header{1,iMeas}.noise.Sli(i),Header{1,iMeas}.noise.Ave(i),Header{1,iMeas}.noise.Rep(i))=2.5*Header{1,iMeas}.noise.timestamp(i);
        end
        Times=squeeze(Times);
        
    else
        Noise=[];
    end
    
end
% cut off spokes for kspace
KSpace = KSpace(:,PAR.dropSpokesTill+1:end,:);
% Building trajectory for radial acquisition from spoke 1 to radial views (= size(kspace,2)/2)
traj = buildRadTraj2D_radial_PC(size(KSpace,1),size(KSpace,2)/2,0,(PAR.dropSpokesTill+2)/2,1,1);
% converting traj to match kspace
traj = cat(1,traj,traj);
traj = reshape(traj,[size(KSpace,1) size(KSpace,2)]);

% displays load
if PAR.verbose
    fprintf('%s from %s has been loaded.\n\n', fname, pathname )
end

FETAL_LOGF (PAR.logf, PAR.verbose, '--- READ RAW DATA MODULE --- end.\n \n')

FETAL_LOGF (PAR.logf, PAR.verbose, '--- RECON PARAM SUMMARY --- start.\n')
params= regexprep(fileread('SET_RECON_PARAMS.m'), '^[%fe].*$', '', 'lineanchors', 'dotexceptnewline');
FETAL_LOGF (PAR.logf, 1, '%s', regexprep(params, '\n\n+', '\n'))
FETAL_LOGF (PAR.logf, PAR.verbose, '--- RECON PARAM SUMMARY --- end.\n \n')

end

