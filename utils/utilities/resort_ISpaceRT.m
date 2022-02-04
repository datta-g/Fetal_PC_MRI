% [IMG,CP0] = resort_ISpaceRT(RT,t,RWaveTimes,nFrames)
%
% Sorts data based on a heart rate model
%
% Original code Christopher Roy 2018
% modifications by Datta Goolaub 2020

function [IMG,CP0] = resort_ISpaceRT(RT,t,RWaveTimes,nFrames)

% compute cardiac phase for RT
CP = Calculate_CardiacPhases(t,RWaveTimes);
CP0=CP; CP=[CP-1;CP;CP+1];
iCP=linspace(0,1,nFrames);

% container for resorted CINE
IMG = zeros([size(RT(:,:,1)),nFrames],'single');

% combines real times into a resorted CINE based on detected cardiac phase
for iFrame=1:nFrames
    x=mygaus([1,iCP(iFrame),mean(diff(iCP))],CP);
    x=x(1:size(RT,3))+x(1+size(RT,3):2*size(RT,3))+x(1+2*size(RT,3):3*size(RT,3));
    IMG(:,:,iFrame) = sum(repmat(permute(x,[2,3,1]),[size(RT,1),size(RT,2),1]).*RT,3)./sum(x(:));
end
end