function res = D(image)

%
% res = D(image)
%
% image = a 2D image
%
% This function computes the finite difference transform of the image
%
% Related functions:
%       adjD , invD 
%
%
% (c) Michael Lustig 2005
% mod DS Goolaub

[sx,sy] = size(image);

% Dx and Dy : x y z t v
Dx = image([2:end,end],:,:,:,:) - image;
Dy = image(:,[2:end,end],:,:,:) - image;


%res = [sum(image(:))/sqrt(sx*sy); Dx(:);  Dy(:)]; 
res = cat(6,Dx,Dy);
% 6th dim is tv