function dmtx = noise_decorrelation_mtx(noise_covariance_mtx)
%
%  dmtx = ismrm_calculate_noise_decorrelation_mtx_from_covariance_mtx(noise_covariance_mtx)
%
%  Computes noise decorrelation matrix for prewhitening from the noise
%  covariance matrix.
%
%  INPUT:
%    - noise_covariance_mtx [coil,coil] : noise covariance matrix
%
%  OUTPUT:
%    - dmtx [coil, coil]                : noise decorrelation matrix
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%
try
    dmtx = inv(chol(noise_covariance_mtx, 'lower'));
catch
    dmtx = eye(size(noise_covariance_mtx));
end
end