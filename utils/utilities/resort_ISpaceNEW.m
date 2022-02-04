% [IMG,CP0] = resort_ISpaceNEW(RT,t,RWaveTimes,nFrames)
%
% Sorts data based on a heart rate model
%
% Original code Christopher Roy 2018

function [IMG,CP0] = resort_ISpaceNEW(RT,t,RWaveTimes,nFrames)
CP = Calculate_CardiacPhases(t,RWaveTimes);
CP0=CP;
CP=[CP-1;CP;CP+1];
iCP=linspace(0,1,nFrames);
IMG = zeros([size(RT(:,:,1)),nFrames],'single');
for iFrame=1:nFrames
    x=mygaus([1,iCP(iFrame),mean(diff(iCP))],CP);
    x=x(1:size(RT,3))+x(1+size(RT,3):2*size(RT,3))+x(1+2*size(RT,3):3*size(RT,3));
    IMG(:,:,iFrame) = sum(repmat(permute(x,[2,3,1]),[size(RT,1),size(RT,2),1]).*RT,3)./sum(x(:));
end
end