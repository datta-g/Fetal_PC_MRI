% [g,e] = Numerical_Gradient_nPar_ISPACE_dsg(iPara,offsets,Times,ISpace,nSpokes,Previous_Para,Dtype)
%
% Computes numerical gradient along a direction for fetal heart rate
%
% Used in entropy minimization using a grid search
%
% Original code Christopher Roy 2018

function [g,e] = Numerical_Gradient_nPar_ISPACE_dsg(iPara,offsets,Times,ISpace,nSpokes,Previous_Para,Dtype)
ScanLength=max(Times(:));
e=zeros(length(offsets),1);
for iOffset=1:length(offsets)
    Para=Previous_Para;
    Para(iPara)=Para(iPara)+offsets(iOffset);
    RW = Linear_Transition_Model(Para,ScanLength);
    RW = Constrained_RWs(RW,ScanLength);
    I = resort_ISpace(ISpace,Times,RW,nSpokes);
    e(iOffset,1)=imagemetric(I,Dtype);
end
f=spline(offsets,smooth(e(:,1)));
fp=fnder(f,1);
g=ppval(fp,0);
end
