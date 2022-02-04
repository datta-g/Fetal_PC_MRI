% NAME :
%           Weights_Radial_Data(traj, normfl)
%
% DESCRIPTION:
%           Sets weights for radial k-space
%
% INPUTS:
%           complex double  traj            matrix holding trajectory (kx + 1i*ky) (samples x profiles x temporal frame)
%           logical         normfl          normalisation flag
%
% OUTPUTS:
%           double          weights         density weights (samples x profiles x temporal frame)
%
% PSEUDOCODE:
%           compute angular separation between spokes
%           compute sector area associated with each sample
%
% NOTES:
%           based on previous code from C. Santelli (2012)
%
% CHANGE LOG:
%     VERSION        DATE             AUTHOR                    CHANGE
%       1.0         2020-11-05       Datta Singh Goolaub        Creation



% Reconstruction code for radial PC MRI using Compressed Sensing, Motion
% Correction and Metric Optimised Gating
%
% Datta Singh Goolaub (2020)
% University of Toronto / The Hospital For Sick Children

function weights = Weights_Radial_Data1(traj, normfl)

% Check input
narginchk(1,2)
if nargin<2, normfl = false; end

% Build containers and extract trajector coordinates
traj_radial_pos = abs(traj(:,1));
weights = zeros(size(traj));

% loops over temporal dimension
for temporal_dim = 1:size(traj,3)
    % get angle of spoke and mask out non spoke data from matrix
    phi  = squeeze(angle(traj(1,:,temporal_dim)));
    zero_mask = sum(abs(traj(:,:,temporal_dim)),1);
    phi(zero_mask==0)=[];%removing zeros in container
    phi(phi<0) = pi + phi(phi<0);
    
    % Sort spokes according to their angles on the interval [0,pi],
    % calculate relative angular distances to neighboring spokes, and
    % finally, get corresponding total relative angles.
    [phi, I] = sort(phi);
    dPhi1   = [phi(2:end) (phi(1)+pi)]-phi; % Left relative angular distance
    dPhi2   = circshift(dPhi1,[0 1]); % Right relative angular distance
    dPhi    = 0.5*(dPhi1+dPhi2);
    w = [];
    
    for i=1:length(phi), w = [w, dPhi(i)*traj_radial_pos];    end
    
    if normfl, w = w./max(w(:)); end
    w(:,I)     = w;
    
    % output weights
    weights(:,(zero_mask~=0),temporal_dim) = w;
    
end

end


