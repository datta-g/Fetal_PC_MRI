function R = ROI_Select_Idea(RT,RR,Scanlength,nFrames,Times,type)
switch type
    case 1
        if numel(RR) ~= 1
            RR = mean(diff(RR));
        end
        MOG_RWaveTimes = Linear_Transition_Model(RR,Scanlength);
    case 2
        %         MOG_RWaveTimes is already in good form
        if min(RR) > 0
            RR = [2*RR(1)-RR(2) RR];
            RR = RR + abs(min(RR));
        end
        MOG_RWaveTimes = RR;
        
end
[R] = resort_ISpace(RT,Times,MOG_RWaveTimes,nFrames);
%  [R] = resortRadKspaceRT_v4(RT,Times, MOG_RWaveTimes, nFrames);
% R = R-Correct_Phase_N(R);

end
% function [IMG,CP0] = resort_ISpace(RT,t,RWaveTimes,nFrames)
% CP = Calculate_CardiacPhases(t,RWaveTimes);
% CP0=CP;
% CP=[CP-1;CP;CP+1];RT=cat(3,RT,RT,RT);
% iCP=linspace(0,1,nFrames);
% IMG = zeros([size(RT(:,:,1)),nFrames],'single');
% for iFrame=1:nFrames
%     x=mygaus([1,iCP(iFrame),mean(diff(iCP))],CP);
%     IMG(:,:,iFrame) = sum(repmat(permute(x,[2,3,1]),[size(RT,1),size(RT,2),1]).*RT,3)./sum(x(:));
% end
% end