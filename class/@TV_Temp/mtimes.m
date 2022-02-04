function res = mtimes(a,b)
% b is x y z t v

if numel(size(b)) == 3
    if a.adjoint
        res = b - circshift(b,[0 0 1]);
    else
        res = b - circshift(b,[0 0 -1]);
    end
else
    if a.adjoint
        res = b - circshift(b,[0 0 0 1 0]);
    else
        res = b - circshift(b,[0 0 0 -1 0]);
    end
end