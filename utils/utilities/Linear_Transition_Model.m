% RWaveTimes =  Linear_Transition_Model(RRs,ScanLength)
%
% Develops heart rate model based on scanlength and RR intervals
%
% Original code Christopher Roy 2018


function RWaveTimes =  Linear_Transition_Model(RRs,ScanLength)
% RWaveTimes is output heart rate model
% RRs is RR interval time
% Scanlength is duration of data acquisition

if length(RRs)==1
    RRs=[RRs,RRs];
end
NP=length(RRs);

if length(RRs)==2 % one param model
    nP1=(ScanLength/NP)/RRs(1);
    nP2=(ScanLength/NP)/RRs(2);
    RRIntervals(1:floor(nP1),1)=RRs(1);
    P1_Diff=(ScanLength/NP)-floor(nP1)*RRs(1);
    Transition=P1_Diff+(1-(P1_Diff/RRs(1)))*RRs(2);
    RRIntervals(end+1,1)=Transition;
    RRIntervals((end+1):(end+1+floor(nP2)),1)=RRs(2);
    RWaveTimes=zeros(length(RRIntervals)+1,1,'single');
    RWaveTimes(2:end)=cumsum(RRIntervals);
    if max(RWaveTimes)<=ScanLength
        RWaveTimes(end+1)=RWaveTimes(end)+(RWaveTimes(end)-RWaveTimes(end-1));
    end
    
elseif length(RRs)>2 % multi param model
    nP1=(ScanLength/NP)/RRs(1);
    RRIntervals(1:floor(nP1),1)=RRs(1);
    P1_Diff=(ScanLength/NP)-floor(nP1)*RRs(1);
    Transition=P1_Diff+(1-(P1_Diff/RRs(1)))*RRs(2);
    RRIntervals(end+1,1)=Transition;
    nP=zeros(NP,1);
    P_Diff=zeros(NP-2,1);
    for p=2:(NP-1)
        nP(p)=((p*ScanLength/NP)-sum(RRIntervals))/RRs(p);
        RRIntervals((end+1):(end+floor(nP(p))),1)=RRs(p);
        P_Diff(p)=(p*ScanLength/NP)-sum(RRIntervals);
        Transition=P_Diff(p)+(1-(P_Diff(p)/RRs(p)))*RRs(p+1);
        RRIntervals(end+1,1)=Transition;
    end
    nP_end=(ScanLength-sum(RRIntervals))/RRs(NP);
    RRIntervals((end+1):(end+ceil(nP_end)),1)=RRs(NP);
    RWaveTimes=zeros(length(RRIntervals)+1,1,'single');
    RWaveTimes(2:end)=cumsum(RRIntervals);
    if max(RWaveTimes)<=ScanLength
        RWaveTimes(end+1)=RWaveTimes(end)+(RWaveTimes(end)-RWaveTimes(end-1));
    end
end

i=find(RWaveTimes>ScanLength,1,'first');
if i<length(RWaveTimes)
    RWaveTimes=RWaveTimes(1:i);
end

if RWaveTimes(end)<=ScanLength
    RWaveTimes(end)=ScanLength+5;
end
end