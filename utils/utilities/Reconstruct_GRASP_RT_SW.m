function [kdatau,ku,wu,indX] = Reconstruct_GRASP_RT_SW(KSpace,ks,ws,nSpokes,wSpokes)
% nSpokes per phase
% wSpokes per overlap between phases
ku=zeros(size(ks,1),nSpokes,'single');
wu=zeros(size(ks,1),nSpokes,'single');
kdatau=zeros(size(KSpace,1),nSpokes,size(KSpace,3),1,'single');
ii=0;
timout = [];
while(1)
    indx = zeros(1,size(KSpace,2));
    ii=ii+1;
    sw=1+(ii-1)*(nSpokes-wSpokes):(ii-1)*(nSpokes-wSpokes)+nSpokes;
    if max(sw)>size(KSpace,2)
        break;
    end
    ku(:,:,ii)=ks(:,sw);
    wu(:,:,ii)=ws(:,sw);
    kdatau(:,:,:,ii)=KSpace(:,sw,:);
    indx(sw)=ii;
    indX(ii,:) = indx;
end
end