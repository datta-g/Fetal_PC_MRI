function res = mtimes(a,bb)

if a(1).op.sensChn==0 || a(1).imageDim(3)>1 %activated if 3D
    
    %     if (a(1).adjoint)
    %         res = gpuNUFFT_adj(a(1).op,bb);
    %     else
    %         res = gpuNUFFT_forw(a(1).op,bb);
    %     end
    %
    %     if (a(1).adjoint)
    %         res=zeros([a(1).imageDim(1:2),size(bb,3),size(bb,4)],'single');
    %     else
    %         res=zeros(a(1).op.params.trajectory_length,a(1).op.sensChn,size(bb,3),size(bb,4),'single');
    %     end
    
    % i is temporal dim
    % j is vel dim
    
    if (a(1).adjoint)
        for i=1:size(bb,3)
            for j=1:size(bb,4)
                count = (i-1)*size(bb,4) + j;
                res(:,:,:,i,j) = gpuNUFFT_adj(a(count).op,bb(:,:,i,j));
            end
        end
    else
        for i=1:size(bb,4)
            for j=1:size(bb,5)
                count = (i-1)*size(bb,4) + j;
                res(:,:,i,j) = gpuNUFFT_forw(a(count).op,bb(:,:,:,i,j));
            end
        end
    end
    
    
    
    %creating data in form of M.*exp(iVo + jVx + ...)
    if (a(1).adjoint)
        res = repmat(abs(mean(res,5)),[1 1 1 1 size(bb,4)]).*exp(1i*angle(res));
    end
    
    
    
    
else
    
    if (a(1).adjoint)
        res=zeros([a(1).imageDim(1:2),size(bb,3),size(bb,4)],'single');
    else
        res=zeros(a(1).op.params.trajectory_length,a(1).op.sensChn,1,size(bb,4),size(bb,5),'single');
    end
    
    for i=1:size(bb,4)      % i is temporal dim
        for j=1:size(bb,5)  % j is vel dim
            count = (i-1)*size(bb,5) + j;
            if (a(count).adjoint)
                res(:,:,:,i,j) = gpuNUFFT_adj(a(count).op,bb(:,:,1,i,j));
            else
                res(:,:,1,i,j) = gpuNUFFT_forw(a(count).op,bb(:,:,:,i,j));
            end
        end
    end
    
    %creating data in form of M.*exp(iVo + jVx + ...)
    if (a(1).adjoint)
%         res = repmat(abs(mean(res,4)),[1 1 1 size(bb,4)]).*exp(1i*angle(res));
        res = repmat(abs(mean(res,5)),[1 1 1 1 size(bb,5)]).*exp(1i*angle(res));

    end
    
    
end

end