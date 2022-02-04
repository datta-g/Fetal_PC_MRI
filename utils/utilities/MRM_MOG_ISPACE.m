% [MOG_RWaveTimes,IMOG,log,offsets,M]=MRM_MOG_ISPACE(ISpace,Times,RRs,rSpokes,Dtype)
%
% Solving for fetal heart rate
%
% Entropy minimization using a grid search
%
% Original code Christopher Roy 2018
% Modifications by Datta Goolaub 2020


function [MOG_RWaveTimes,IMOG,log,offsets,M]=MRM_MOG_ISPACE(ISpace,Times,PAR,RT_del)

Scanlength=range(Times(:));
% PAR.MOG.RRrange = PAR.MOG.RRrange(70:112);
% solve for initial single parameter model
M=zeros(length(PAR.MOG.RRrange),1);
for iRR=1:length(PAR.MOG.RRrange)
    I = resort_ISpaceRT(ISpace,Times,Linear_Transition_Model(PAR.MOG.RRrange(iRR),Scanlength),PAR.MOG.CardPhase);
    M(iRR,1)=imagemetric(I,{PAR.MOG.METRIC});
    S(:,:,:,1,iRR)=resort_ISpaceRT(RT_del(:,:,:,1),Times,Linear_Transition_Model(PAR.MOG.RRrange(iRR),Scanlength),PAR.MOG.CardPhase);%del
    S(:,:,:,2,iRR)=resort_ISpaceRT(RT_del(:,:,:,2),Times,Linear_Transition_Model(PAR.MOG.RRrange(iRR),Scanlength),PAR.MOG.CardPhase);%del

end

% finds average heart rate with min entropy
PreviousPara = PAR.MOG.RRrange(find(M==min(M),1,'first'));
% develops single parameter model for heart rate
MOG_RWaveTimes = Constrained_RWs(Linear_Transition_Model(PreviousPara,Scanlength),Scanlength);
offsets=[];
log=struct('RWaveTimes',[],'Entropy',[],'ES',[]);
log(1).RWaveTimes=MOG_RWaveTimes;
log(1).Entropy=min(M);
log(1).ES=M;

% number of terms in model
NP=floor(log2(length(MOG_RWaveTimes)));

% solve for multiparameter model by doubling parameters in each stage
for loop=1:NP
    % solve for new heart rate model with additional parameters
    PreviousPara = Double_Parameters(PreviousPara);
    [~,PreviousPara,ES,offsets] = Gradient_MOG_nPara_ISPACE(Times,ISpace,PAR.MOG.CardPhase,PreviousPara,{PAR.MOG.METRIC});
    MOG_RWaveTimes = Constrained_RWs(Linear_Transition_Model(PreviousPara,Scanlength),Scanlength);

    % resorts image for logging
    I = resort_ISpaceRT(ISpace,Times,MOG_RWaveTimes,PAR.MOG.CardPhase);
    log(loop+1).RWaveTimes=MOG_RWaveTimes;
    log(loop+1).Entropy=imagemetric(I,{PAR.MOG.METRIC});
    log(loop+1).ES=ES;
    log(loop+1).reorderedCINE = I;
end
% output heart rate model
MOG_RWaveTimes = Constrained_RWs(Linear_Transition_Model(PreviousPara,Scanlength),Scanlength);

if PAR.verbose
    display('Final Reconstruction')
end
IMOG = resort_ISpaceRT(ISpace,Times,MOG_RWaveTimes,PAR.MOG.CardPhase);
MOG_RWaveTimes = MOG_RWaveTimes';

end

