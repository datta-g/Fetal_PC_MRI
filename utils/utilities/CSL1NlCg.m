% res = CSL1NlCg(param)
%
% Compressed sensing reconstruction of undersampled k-space MRI data
%
% L1-norm minimization using non linear conjugate gradient iterations
%
% Given the acquisition model y = E*x, and the sparsifying transform W,
% the pogram finds the x that minimizes the following objective function:
%
% f(x) = ||E*x - y||^2 + lambda * ||W*x||_1
%
% Based on the paper: Sparse MRI: The application of compressed sensing for rapid MR imaging.
% Lustig M, Donoho D, Pauly JM. Magn Reson Med. 2007 Dec;58(6):1182-95.
%
% Ricardo Otazo, NYU 2008

% modifications by Datta Goolaub 2020

function [x] = CSL1NlCg(rec0,param)

% starting point
x=rec0;

% line search parameters
maxlsiter = 150 ;
gradToll = 1e-4 ;
param.l1Smooth = 1e-15;
alpha = 0.01;
beta = 0.6;
t0 = 1;
ite_count = 0;

% compute g0  = grad(f(x))
g0 = grad(x,param);
dx = -g0;

% iterate
while(1)
    reset(parallel.gpu.GPUDevice.current())
    % preobj
    preF = preobjective(x,dx,param);
    % backtracking line-search
    f0 = objective(preF,0,param);
    t = t0;
    
    f1 = objective(preF,t,param);
    
    lsiter = 0;
    while (f1.Obj > f0.Obj - alpha*t*abs(g0(:)'*dx(:)))^2 & ...
            (lsiter<maxlsiter)
        lsiter = lsiter + 1;
        t = t * beta;
        f1 = objective(preF,t,param);
    end
    
    
    if lsiter == maxlsiter
        disp('Error - line search ...');
        return;
    end
    % control the number of line searches by adapting the initial step search
    if lsiter > 2
        t0 = t0 * beta;
    end
    if lsiter<1
        t0 = t0 / beta;
    end
    
    % update x
    x = (x + t*dx);
    
    if param.verbose
        c_snr = abs(x(round(0.45*size(x,1)):round(0.55*size(x,1)),round(0.45*size(x,2)):round(0.55*size(x,2)),round(size(squeeze(x),3))));
        c_snr = mean(c_snr(:))/std(c_snr(:));
        fprintf('ite = %d, cost = %f, csnr = %f \n',ite_count,f1.Obj, c_snr);
    end
    
    %conjugate gradient calculation
    g1 = grad(x,param);
    bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
    
    g0 = g1;
    dx =  - g1 + bk*dx;
    
    ite_count = ite_count + 1;
    
    if (ite_count > param.nite) || (norm(dx(:)) < gradToll)
        break;
    end
    
end
return;
end

function res = preobjective(x,dx,param) %**********************************
l2x = (param.GPU*x).*param.oneM;
l2dx = (param.GPU*dx).*param.oneM;

res.L2x = l2x; res.L2dx = l2dx;

Objx=zeros(size(x,1),size(x,2),size(x,3),size(x,4),size(x,5),2,length(param.Weights));
Objdx=Objx;
if param.Weights(1)
    for i = 1:length(param.Weights)
        s= param.Reg{i}*x;
        Objx(:,:,:,:,:,1:size(s,6),i) = s;
        Objdx(:,:,:,:,:,1:size(s,6),i) = param.Reg{i}*dx;
    end
end
res.x = Objx; res.dx = Objdx;
end



function res = objective(pres, t, param) %**********************************

% L2-norm part
wx = pres.L2x + t*pres.L2dx - param.y;
L2Objx=sum(abs(wx(:)).^2);
res.reg=[];
SumRegs = 0;
if param.Weights(1)
    for i = 1:length(param.Weights)
        w = (pres.x(:,:,:,:,:,:,i) + t*pres.dx(:,:,:,:,:,:,i));
        res.reg(i) = sum((w(:).*conj(w(:))+param.l1Smooth).^(1/2));
        SumRegs = SumRegs + param.Weights(i)*res.reg(i);
    end
end
res.reg(end+1) = L2Objx;
res.Obj=L2Objx + SumRegs;

end

function g = grad(x,param)%***********************************************

% L2-norm part
L2Gradx = 2*(param.GPU'*((param.GPU*x).*param.oneM -param.y));

Grad = zeros(size(x,1),size(x,2),size(x,3),size(x,4),size(x,5));
if param.Weights(1)
    for i = 1:length(param.Weights)
        wx = param.Reg{i}*x;
        Grad = Grad + param.Weights(i)*(param.Reg{i}'*(wx.*(wx.*conj(wx)+param.l1Smooth).^(-0.5)));
    end
end

% complete gradient
g=L2Gradx+Grad;
end