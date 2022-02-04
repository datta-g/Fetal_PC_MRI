%  RW = Constrained_RWs(RW,ScanLength)
%
% creates a heart rate model
%
% original code Christopher Roy 2018

function RW = Constrained_RWs(RW,ScanLength)
% RW is constrained heart rate model
% RW is input heart rate model
% ScanLength is the data acquisition time

RW = makerow(RW);

if RW(1)~=0
    RW=[0,RW];
end


if mean(diff(RW))<750
    i=find(abs(diff(RW,2))>35);
elseif mean(diff(RW))<1500
    i=find(abs(diff(RW,2))>150);
else
    i=find(abs(diff(RW,2))>500);
end

if ~isempty(i);
    RR=diff(RW);
    RR(i+1)=RR(i+1)-0.5*(RR(i+1)-RR(i));
    RR(i)=RR(i)+0.5*(RR(i+1)-RR(i));
    RW=cumsum([0,RR]);
end


while RW(end)<ScanLength
    RW(end)=RW(end)+ScanLength+10;
end


end

