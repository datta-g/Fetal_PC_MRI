function [Ksp] = Translate_Correct_KSpace(Ksp1,offsets,k)
Ksp = zeros(size(Ksp1));

for vel = 1:size(Ksp1,4) 
    for chn = 1:size(Ksp1,3)
        for ech = 1:size(Ksp1,2)
            Ksp(:,ech,chn,vel) = Ksp1(:,ech,chn,vel).*exp(2*pi*1i*(offsets(1,ech)*real(k(:,ech))+offsets(2,ech)*imag(k(:,ech))));
        end
    end
end

end