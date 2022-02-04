function I = resort_ISpace(ISpace,Times,RW,rSpokes)
sTim = Times;
% CP = Calculate_CardiacPhases(Times,RW);I = Cardiac_Bins_Sliding_Window_ISpace(CP,ISpace,size(ISpace,3)/rSpokes,0.5);

% 
I = resort_ISpaceNEW(ISpace,Times,RW,rSpokes);

% [ISpace, sTim] = interpolateAlong3rdDim(ISpace,Times,2.7);
% [ISpace, sTim] = interpolateAlong3rdDim(ISpace,Times,34);
% [I] = resortRadKspaceRT_v4(ISpace,sTim, RW, rSpokes);

 
 
% CP = Calculate_CardiacPhases(Times,RW);
% iCP=linspace(0,1,ceil(size(ISpace,3)/rSpokes));
% 
% CP=cat(1,CP-1,CP,CP+1);
% I=Interp_Phantom(cat(3,ISpace,ISpace,ISpace),CP,iCP);



% I = resort_ISpace_SW(ISpace,k,Times,RW,rSpokes,0);

% Times=Times-min(Times(:));
% rFrames=floor(size(ISpace,3)/rSpokes);
% if ~isempty(RW)
% rPhases = Calculate_CardiacPhases(Times,RW);
% [~,i]=sort(rPhases(:));
% else
% [~,i]=sort(Times(:));
% end
% I=zeros(size(ISpace,1),size(ISpace,2),rFrames,'single');
% for iFrame=1:rFrames
% D = local_angular_sampling_density(k(:,i((iFrame-1)*rSpokes+1:iFrame*rSpokes)));
% I(:,:,iFrame)=sum(ISpace(:,:,i((iFrame-1)*rSpokes+1:iFrame*rSpokes)).*repmat(D,[size(ISpace,1),size(ISpace,2),1]),3);
% end

end



% function I = resort_ISpace(ISpace,k,Times,RW,rFrames)
% if length(Times)>600
%     nSpokes=400;
% elseif length(Times)>300
%     nSpokes=200;
% else
%     nSpokes=100;
% end
% rSpokes=floor(size(ISpace,3)/rFrames);
% rPhases = Calculate_CardiacPhases(Times,RW);
% [~,i]=sort(rPhases(:));
% I=zeros(128,128,rFrames,'single');
% for iFrame=1:rFrames
%     w=1+((iFrame-1)*rSpokes):((iFrame-1)*rSpokes)+nSpokes;
%     w(w>size(ISpace,3))=w(w>size(ISpace,3))-size(ISpace,3);
%     D = local_angular_sampling_density(k(:,i(w)));
%     I(:,:,iFrame)=sum(ISpace(:,:,i(w)).*repmat(D,[size(ISpace,1),size(ISpace,1),1]),3);
% end
% end


