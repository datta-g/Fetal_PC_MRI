%=========================================================================
% Create radial trajectory for NUFFT
function [k,ANGLES] = buildRadTraj2D_radial_PC(no_samples, no_profiles, alt_prof, STARTOFFSET, gafl, normfl, dim_t, t_offset, dim_z, z_offset)
%==========================================================================
% Returns radial trajectory in a complex valued data array. The real part
% corresponds to the x- and the imaginary part to y-component respectively.
%
% Inputs:
% -------
% no_samples:   Number of samples along each projection.
% no_rofiles:   Number of projections per frame and slice.
% alt_prof:     Alternating profiles flag.
% gafl:         Golden angle flag.
% normfl:       Coordinates normalization flag.
% dim_t:        Number of time frames.
% t_offset:     Profile offset between two adjacent time frames.
% dim_z:        Number of slices in z-direction.
% z_offset:     Profile offset between two adjacent slices.
%
% Outputs:
% --------
% k:            k-space coordinates.
%
% Function calls:    none
% ---------------
%
% Claudio Santelli, 11/06/2012
%==========================================================================

% Check input
%--------------------------------------------------------------------------
narginchk(3,9)
if nargin<10 || isempty(z_offset), z_offset = 0; end
if nargin<9 || isempty(dim_z),    dim_z    = 1; end
if nargin<8 || isempty(t_offset), t_offset = 0; end
if nargin<7 || isempty(dim_t),    dim_t    = 1; end
if nargin<6 || isempty(normfl),   normfl   = true; end
if nargin<5 || isempty(gafl),     gafl     = false; end

% Build trajectory
%--------------------------------------------------------------------------
k = [];

% Initial spoke along ky-axis
k0 = [zeros(1,no_samples); linspace(-(no_samples/2),(no_samples/2),no_samples)];

% Angle increment
if gafl
    goldenRatio = (sqrt(5)+1)/2;
    dPhi        = pi/goldenRatio;
%     goldenRatio = (3-sqrt(5));
%     dPhi        = pi*goldenRatio;
else
    dPhi = pi/no_profiles;
end
for z=1:dim_z
    for t=1:dim_t
        for i=1:no_profiles
            % Update rotation matrix
            rot_angle = (STARTOFFSET-1)*dPhi + ((i-1)+(t-1)*t_offset+(z-1)*z_offset)*dPhi;
            if alt_prof && ~mod(i,2)
                rot_angle = rot_angle+pi;
            end
            ANGLES(i,z,t)=mod(rot_angle,pi);
            R = [cos(rot_angle), -sin(rot_angle);
                sin(rot_angle),  cos(rot_angle)];
            % Rotate k0 vector accordingly
            ktmp       = (R*k0).';
            k(:,i,z,t) = ktmp(:,1)+1i*ktmp(:,2);
        end
    end
end

if normfl, k = k./no_samples; end

end


