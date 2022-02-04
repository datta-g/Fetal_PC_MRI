
function f = mygaus(param,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function used with leastsq to fit data to
% SF=height * exp(-(x-centre)^2/(2o^2)))
% param(1) = height
% param(2) = center
% param(3) = sigma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f= param(1)*exp((-(x-param(2)).^2)/(2*param(3)^2));
end

function I = Decimate_CINE_v2(I,n)
if n>1
    for t=1:size(I,3)
        for u=1:size(I,4)
            I(:,:,t,u) = Decimate_Resolution_v2(I(:,:,t,u),n);
        end
    end
end
end

function MDR = Decimate_Resolution_v2(z,n)
z=double(abs(z));
kern = ones(n,n)*1/n;
MDR = conv2(z,kern);
MDR = MDR(n:end, n:end);
end
