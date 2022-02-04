function I = Window_CINE(I,mn,mx)
I=abs(I-min(I(:)));
I=I./abs(max(I(:)));
I(I<mn*max(I(:)))=mn*max(I(:));
% I(I<mn*min(I(:)))=mn*min(I(:));
I(I>mx*max(I(:)))=mx*max(I(:));
I=I./max(abs(I(:)));
end


