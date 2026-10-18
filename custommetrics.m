function M = custommetrics(testI, testM, trueI, trueM, dim)
% CUSTOMMETRICS calculates the Dice coefficient, Intersection over Union,
% and inverted Root Mean Squared Error for the given input data

% Determine dimensionality
if dim == 2
    dims = [1 2];
    denom = size(testI,1)*size(testI,2);
elseif dim == 3
    dims = [1 2 3];
    denom = size(testI,1)*size(testI,2)*size(testI,3);
else
    error('Incorrect dimension number.');
end

% Calculate iRMSE
iRMSE = -sqrt(squeeze(sum((trueI - testI).^2,dims))./denom);

% dice & iou parameters
area_i = squeeze(sum(testM.*trueM, dims));
area_t = (squeeze(sum(testM, dims)) + squeeze(sum(trueM, dims)));

% dice
dsc = 2.*area_i./area_t;

% % iou
% iou = area_i./(area_t - area_i);
% 
% M = [dsc, iou, iRMSE];
M = [dsc, iRMSE];
end