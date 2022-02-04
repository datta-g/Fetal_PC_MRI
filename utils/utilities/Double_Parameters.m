% Para = Double_Parameters(PreviousPara)
% 
% doubles the number of parameters in model
% 
% Original code Christopher Roy 2018

function Para = Double_Parameters(PreviousPara)
Para=[];
for loop=1:length(PreviousPara)
    Para=[Para,PreviousPara(loop),PreviousPara(loop)];
end
end
