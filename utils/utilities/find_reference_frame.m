function [frame_ref, sumOfCorrel,CorMat] = find_reference_frame(time_series, ROI)

time_series = round(double(time_series));
% preparing ROI coordinates struct
switch nargin
    case 2
        % do nothing
    case 1
        ROI.y = 1:size(time_series,1);
        ROI.x = 1:size(time_series,2);
end


% cropping to selected region of non-zeros
time_series=time_series(ROI.y,ROI.x,:);

% creating a correlation matrix between all frames
for cur_frame = 1:size(time_series,3)
    for other_frame = (cur_frame+1):size(time_series,3)
        curentFrame = time_series(:,:,cur_frame);
        otherFrame = time_series(:,:,other_frame);
        CorMat(cur_frame,other_frame) = nmi(curentFrame(:),otherFrame(:));
        CorMat(other_frame,cur_frame) = CorMat(cur_frame,other_frame);
    end
end
sumOfCorrel = sum(CorMat,1);
frame_ref = find(sumOfCorrel == max(sumOfCorrel));