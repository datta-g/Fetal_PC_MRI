function res = ctranspose(a)
for i=1:length(a)
a(i).adjoint = xor(a(i).adjoint,1);
res = a;
end

