function dsc = dicecustom(~, m_CT, ~, m_mri, dim)
% DICECUSTOM calculates the Dice coefficient for the given input masks
%   dsc = DICECUSTOM(~, m_CT, ~, m_mri, dim)

% Determine dimensionality
if dim == 2
    dims = [1 2];
elseif dim == 3
    dims = [1 2 3];
else
    error('Incorrect dimension number.');
end

% dice parameters
area_i = squeeze(sum(m_CT.*m_mri, dims));
area_t = (squeeze(sum(m_CT, dims)) + squeeze(sum(m_mri, dims)));

% dice
dsc = 2.*area_i./area_t;
end