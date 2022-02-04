% [REF,PreviousPara,ES,offsets] = Gradient_MOG_nPara_ISPACE(Times,ISpace,nSpokes,PreviousPara,Dtype)
%
% Solving for fetal heart rate
%
% Entropy minimization using a grid search
%
% Original code Christopher Roy 2018

function [REF,PreviousPara,ES,offsets] = Gradient_MOG_nPara_ISPACE(Times,ISpace,nSpokes,PreviousPara,Dtype)

nPara=length(PreviousPara);
ScanLength=max(Times(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Gradient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
StepSizes=[0 1 2 5 10 15];
G=zeros(nPara,1);

trial=0;
while trial<5
    
    nOffsets=15;
    
    steps=5;offsets=-0.5*(nOffsets-1)*steps:steps:0.5*(nOffsets-1)*steps;clear nOffsets steps;
    ES=zeros(length(offsets),nPara);
    
    for iPara=1:nPara
        [g,e] = Numerical_Gradient_nPar_ISPACE_dsg(iPara,offsets,Times,ISpace,nSpokes,PreviousPara,Dtype);
        G(iPara,1)=g;
        ES(:,iPara)=e;
    end
    
    [~,i]=min(ES);
    if isempty(find(offsets(i)~=0, 1))
        break;
    end
    
    
    %%% Minimum Jitter
    [i,j]=find(ES==min(ES(:)),1,'first');
    Para=PreviousPara;
    Para(j)=Para(j)+offsets(i);
    RW = Linear_Transition_Model(Para,ScanLength);
    RW = Constrained_RWs(RW,ScanLength);
    ParaJitter=Para;
    img = resort_ISpace(ISpace,Times,RW,nSpokes);
    MinimumJitter=imagemetric(img,Dtype);
    
    %%% Minimum Global Shift
    [~,i]=min(ES);
    Para=PreviousPara;
    Para=makerow(Para)+offsets(i);
    RW = Linear_Transition_Model(Para,ScanLength);
    RW = Constrained_RWs(RW,ScanLength);
    ParaGlobal=Para;
    img = resort_ISpace(ISpace,Times,RW,nSpokes);
    MinimumGlobal=imagemetric(img,Dtype);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Line SearchovSpks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    LineEntropy=zeros(length(StepSizes),1);
    LinePara=zeros(nPara,length(StepSizes));
    
    parfor loop=1:length(StepSizes)
        Para=PreviousPara;
        Shifts=floor(StepSizes(loop)*G./max(abs(G)));
        Para=makerow(Para)-makerow(Shifts);
        RW = Linear_Transition_Model(Para,ScanLength);
        if loop>1
            RW = Constrained_RWs(RW,ScanLength);
        end
        LinePara(:,loop)=Para;
        img = resort_ISpace(ISpace,Times,RW,nSpokes);
        tempE=imagemetric(img,Dtype);
        LineEntropy(loop,1)=tempE;
    end
    
    [MinimumGrad,i]=min(LineEntropy);
    ParaGrad=LinePara(:,i);
    
    [~,i]=min([MinimumJitter,MinimumGlobal,MinimumGrad]);
    if i==1
        PreviousPara=ParaJitter;
    elseif i==2
        PreviousPara=ParaGlobal;
    elseif i==3
        PreviousPara=ParaGrad;
    end
    
    trial=trial+1;
end
MOG_RWaveTimes = Constrained_RWs(Linear_Transition_Model(PreviousPara,ScanLength),ScanLength);
REF = resort_ISpace(ISpace,Times,MOG_RWaveTimes,nSpokes);
end
