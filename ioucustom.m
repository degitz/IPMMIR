function iou = ioucustom(~, m_CT, ~, m_mri, dim)
% IOUCUSTOM calculates the intersection over union for the given input masks
%   iou = IOUCUSTOM(~, m_CT, ~, m_mri, dim)

% Determine dimensionality
if dim == 2
    dims = [1 2];
elseif dim == 3
    dims = [1 2 3];
else
    error('Incorrect dimension number.');
end

% iou parameters
area_i = squeeze(sum(m_CT.*m_mri, dims));
area_t = (squeeze(sum(m_CT, dims)) + squeeze(sum(m_mri, dims)));

% iou
iou = area_i./(area_t - area_i);
end