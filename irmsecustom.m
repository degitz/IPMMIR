% Calculate inverted root mean squared error
function iRMSE = irmsecustom(i_CT, ~, i_mri, ~, dim)
% IRMSECUSTOM calculates the inverted root mean squared error for the given image slices.
%   iRMSE = IRMSECUSTOM(i_CT, ~, i_mri, ~, dim)

% Determine dimensionality
if dim == 2
    dims = [1 2];
    iRMSE = -sqrt(squeeze(sum((i_mri - i_CT).^2,dims))./(size(i_CT,1)*size(i_CT,2)));
elseif dim == 3
    dims = [1 2 3];
    iRMSE = -sqrt(squeeze(sum((i_mri - i_CT).^2,dims))./(size(i_CT,1)*size(i_CT,2)*size(i_CT,3)));
else
    error('Incorrect dimension number.');
end


end