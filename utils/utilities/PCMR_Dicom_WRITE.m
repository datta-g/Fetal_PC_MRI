% Reconstruction code for radial PC MRI using Compressed Sensing, Motion
% Correction and Metric Optimised Gating
%
% Datta Singh Goolaub (2018)
% University of Toronto / The Hospital For Sick Children

function PCMR_Dicom_WRITE(writePath,subname, fname,pcmri_recon,Header, SeriesNumber, RRs, resftraction)

% % check if path exists
% if ~exist([writePath,'/DicomsA'],'dir')
%     mkdir([writePath,'/DicomsA'])
% end

% resize image to correct for changed res
for i = 1:size(pcmri_recon,4)
    for j = 1:size(pcmri_recon,5)
        pcmri_recon1(:,:,:,i,j) = imresize(pcmri_recon(:,:,:,i,j),1/resftraction);
    end
end
pcmri_recon = pcmri_recon1; clear pcmri_recon1;

% Patient=Header{end}.hdr.Dicom.tPatientName;
Patient.FamilyName = 'fetus Sheep';
Patient.GivenName = subname;
Date = subsref(strsplit(Header{end}.hdr.Config.ExamMemoryUID,'_'),struct('type','()','subs',{{4}}));
Date =  Date{1}(1:8);
info.StudyDate = Date;
info.SeriesDate = Date;
info.AcquisitionDate = Date;
info.ImageDate = Date;
info.InstitutionName = Header{end}.hdr.Dicom.InstitutionName;
info.StudyInstanceUID=Header{end}.hdr.Config.FrameOfReference;
info.SeriesInstanceUID=dicomuid;
info.FrameOfReferenceUID = Header{end}.hdr.Config.FrameOfReference;
info.Manufacturer='SIEMENS';
info.PixelSpacing= repmat(Header{end}.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV/Header{end}.hdr.Dicom.lBaseResolution,[1 2]);
info.Modality = 'MR';
info.PatientName=Patient;
info.ImageComments='N/A';
info.PatientComments='N/A';
SlicePos=Header{end}.image.slicePos(4:7,1);
info.SliceLocation=Header{end}.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dThickness;
info.StudyDate=Date;
info.SeriesDate=Date;
info.AcquisitionDate=Date;
info.ContentDate=Date;
venc = Header{end}.hdr.MeasYaps.sAngio.sFlowArray.asElm{1};
info.ImagePositionPatient=SlicePos(1:3,1);
info.Private_0019_1015=SlicePos(1:3,1);
info.Private_0019_1016=SlicePos(1:3,1);

slicePos = SlicePos;
a = slicePos(1); b = slicePos(2); c = slicePos(3); d = slicePos(4);
R = [a*a+b*b-c*c-d*d, 2*b*c-2*a*d, 2*b*d+2*a*c;
    2*b*c+2*a*d, a*a-b*b+c*c-d*d, 2*c*d-2*a*b;
    2*b*d-2*a*c, 2*c*d+2*a*b, a*a-b*b-c*c+d*d  ];
newR = R*[1 0 0; 0 1 0]';
info.ImageOrientationPatient = newR(:);
info.CardiacNumberOfImages = size(pcmri_recon,3);
nominalInterval = round(mean(RRs));
info.NominalInterval = nominalInterval;
frame_array = 0:nominalInterval/size(pcmri_recon,3):nominalInterval-nominalInterval/size(pcmri_recon,3);

protocolName = Header{end}.hdr.Dicom.tProtocolName;
outPathM = [writePath,'/Dicoms/' fname '_M/']; mkdir(outPathM);
outPathP = [writePath,'/Dicoms/' fname '_P/']; mkdir(outPathP);
iCount = 0;
data_I = pcmri_recon(:,:,:,1).*conj(pcmri_recon(:,:,:,2));
mag = sqrt(abs(data_I));
mag = uint16(mag/max(mag(:))*4096);
flow = angle(data_I);
flow = uint16(flow/pi*2048 + 2048);
for imType = 1:2
    if imType == 1 % fill magnitude or venc header info
        I = mag;
        info.SequenceName = 'fl2d1r2';
        info.ImageType = 'ORIGINAL\PRIMARY\M\RETRO\DIS2D';
        info.SeriesNumber=SeriesNumber;
        info.Private_0051_0016 = 'p2 M/RETRO/DIS2D';
        outPath = outPathM;
    else
        I = flow;
        info.SequenceName = ['fl2d1_v' num2str(venc.nVelocity) 'in'];
        info.ImageType = 'DERIVED\PRIMARY\P\RETRO\DIS2D';
        info.SeriesNumber=SeriesNumber+1;
        info.Private_0051_0014 = ['v' num2str(venc.nVelocity) '_through'];
        info.Private_0051_0016 = 'p2 P/RETRO/DIS2D';
        outPath = outPathP;
    end
    for iDCM=1:size(pcmri_recon,3)
        FrameCounter=iDCM;
        info.TriggerTime = frame_array(iDCM);
        
        iCount = iCount + 1;
        info.InstanceNumber=iCount;
        
        info.SliceThickness=Header{end}.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dThickness;
        if imType ==1
            info.Filename=[outPath, protocolName, '_Magnitude_Frame_',num2str(FrameCounter),'.dcm'];
        else
            info.Filename=[outPath, protocolName, '_Phase_Frame_',num2str(FrameCounter),'.dcm'];
        end
        dicomwrite(I(:,:,FrameCounter),info.Filename,info);
    end
    iCount = 0;
end
clearvars -except FileList Path Data
clc
display('___________________________________________________________________________')