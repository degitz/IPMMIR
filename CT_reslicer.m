function [i_CT, m_CT, varargout] = CT_reslicer(i_ct, m_ct, r, inp4, inp5)
% ct - ct data
% mri - mri data
% r - number of ct slices per mri slice (mri res / ct res)

% Determine indices of start/stop for mri slices in ct
ibase = 0:r:size(i_ct,3); % base starts & stops
ifull(1,:) = ceil(ibase(1:end-1))+1; % first full slices
ifull(2,:) = floor(ibase(2:end)); % last full slices

% Find non-integer amounts contributing to slices
ipart(1,:) = ceil(ibase(1:end-1))-ibase(1:end-1);
ipart(2,:) = ibase(2:end)-floor(ibase(2:end));

% Determine what to do with input
if size(inp4, 1) == 1
    Bpix = inp4(1); Bslc = inp4(2); 
    
    % Isolate first & last slice
    firstslice = 1 + Bslc;
    lastslice = inp5;

    % Reslice data
    i_CT = SliceCombine(i_ct, ifull+Bpix, ipart);
    m_CT = SliceCombine(m_ct, ifull+Bpix, ipart);
    i_CT = i_CT(:,:,firstslice:lastslice);
    m_CT = m_CT(:,:,firstslice:lastslice);
else
    i_mri = inp4; m_mri = inp5;

    % Determine size difference of data
    slcdiff = size(ifull,2)-size(i_mri,3);
    pixdiff = size(i_ct,3)-ifull(2,end);

    % Preallocate data
    MTpix = zeros(pixdiff+1, 2);
    Bslc = zeros(pixdiff+1, 1);

    % Compare CT slice assignments with available subslice offsets. Last offset
    % is excluded, since partial pixel addition would extend past matrix bounds.
    for pixofs = 0:pixdiff
        % fprintf('\nTesting a subslice offset of %i pixel(s).\n', pix)
        i_CT = SliceCombine(i_ct, ifull+pixofs, ipart);
        m_CT = SliceCombine(m_ct, ifull+pixofs, ipart);

        % For extra slices in CT data, determine offset with MRI
        % Preallocate data
        MTslc = zeros(slcdiff+1, 2);
        for slcofs = 0:slcdiff

            % Isolate current number of slices
            i_CT_ofs = i_CT(:,:,1+slcofs:end-slcdiff+slcofs);
            m_CT_ofs = m_CT(:,:,1+slcofs:end-slcdiff+slcofs);

            MTslc(1+slcofs,:) = custommetrics(i_CT_ofs, m_CT_ofs, i_mri, m_mri, 3);
        end

        % Mostly vectorized version of above code. It is less efficient.
        % % For extra slices in CT data, determine offset with MRI
        % % Preallocate data
        % met_ofs = zeros(size(i_CT,3)-size(mri,3)+1,size(i_CT,3)-size(mri,3));
        % CT_ofs = zeros(size(i_CT,1), size(i_CT,2), size(mri,3), size(i_CT,3)-size(mri,3));
        % for slc = 0:size(i_CT,3)-size(mri,3)
        %     % Isolate current number of slices
        %     CT_ofs(:,:,:,slc+1) = i_CT(:,:,1+slc:end-(size(i_CT,3)-size(mri,3))+slc);
        % end
        %
        % % RMSE
        % RMSE = sqrt(squeeze(sum((mri - CT_ofs).^2,[1 2]))./(size(CT_ofs,1)*size(CT_ofs,2)));
        %
        % % dice parameters
        % CTm = CT_ofs > 0;
        % area_i = squeeze(sum(CTm.*m_mri, [1 2]));
        % area_t = (squeeze(sum(CTm, [1 2])) + squeeze(sum(m_mri, [1 2])));
        %
        % % dice
        % dsc = 2.*area_i./area_t;
        %
        % met_ofs = mean(dsc - RMSE);

        % Find best
        Bslc(1+pixofs) = IsolateBest(MTslc)-1;

        % Isolate best from the rest
        i_CT = i_CT(:,:,1+Bslc(1+pixofs):end-slcdiff+Bslc(1+pixofs));
        m_CT = m_CT(:,:,1+Bslc(1+pixofs):end-slcdiff+Bslc(1+pixofs));

        % fprintf('For a %i pixel subslice offset, a %i slice offset was best.\n', pix, Bslc(1+pix))

        % Evaluate
        MTpix(1+pixofs,:) = custommetrics(i_CT, m_CT, i_mri, m_mri, 3);
    end

    % Find best
    Bpix = IsolateBest(MTpix);

    % Isolate best from the rest
    Bslc = Bslc(Bpix); Bpix = Bpix-1;

    % Isolate correct slice & pixel values by adding the best offsets
    firstslice = 1 + Bslc;
    lastslice = size(i_mri,3) + Bslc;
    % firstpixel = ifull(1,firstslice) + Bpix;
    % lastpixel = ifull(2,lastslice) + Bpix;
    % fprintf('\nOverall best results obtained with offsets of %i pixel(s) and %i slice(s).\n', Bpix, Bslc)
    % fprintf('This includes slices %i to %i and pixels %i to %i, plus %.4f and %.4f of pixels %i and %i, respectively.\n', firstslice, lastslice, firstpixel, lastpixel, ipart(1,firstslice), ipart(2,lastslice), firstpixel-1, lastpixel+1)
    i_CT = SliceCombine(i_ct, ifull+Bpix, ipart);
    m_CT = SliceCombine(m_ct, ifull+Bpix, ipart);
    i_CT = i_CT(:,:,firstslice:lastslice);
    m_CT = m_CT(:,:,firstslice:lastslice);

    % Set necessary variables as outputs
    varargout{1} = Bpix;
    varargout{2} = Bslc;
end
end

% Determine best result
function B = IsolateBest(met)
nume = 10; % number of extrema before combination

% Find best
[~, met_bi1] = maxk(met(:,1), nume);
[~, met_bi2] = maxk(met(:,2), nume);

% Combine all best into one array
B = intersect(met_bi1, met_bi2, 'stable');

% Check if none meet requirements
if isempty(B)

    % Create array of li values for best mi values
    met_tmp = nan(size(met(:,2))); met_tmp(met_bi1) = met(met_bi1,2);

    % Find max
    [~, B] = max(met_tmp);
end
B = B(1);
end

% Combine microCT slices into MRI slices
function IMG = SliceCombine(img, ifull, ipart)
% Preallocate data
IMG = zeros(size(img,1), size(img,2), size(ifull,2), 'like', img);

% Check if on last pixel slice offset
if ifull(2,end) == size(img,3)
    % Flip partial slices so we go in reverse
    ipart = fliplr(ipart);

    % Update loop parameters
    ystart = 1;
    ystop = size(ifull,2)-1;

    % Hard code last iteration (last slice should end at an integer)
    IMG(:,:,end) = sum(img(:,:,ifull(1,end):ifull(2,end)), 3) + img(:,:,ifull(1,end)-1).*ipart(1,end);
else
    % Update loop parameters
    ystart = 2;
    ystop = size(ifull,2);

    % Hard code first iteration (first slice should start at an integer)
    IMG(:,:,1) = sum(img(:,:,ifull(1,1):ifull(2,1)), 3) + img(:,:,ifull(2,1)+1).*ipart(2,1);
end
for y = ystart:ystop
    % Sum full slices + part of previous slice + part of next slice
    IMG(:,:,y) = sum(img(:,:,ifull(1,y):ifull(2,y)), 3) + img(:,:,ifull(1,y)-1).*ipart(1,y) + img(:,:,ifull(2,y)+1).*ipart(2,y);
end

% Normalize
IMG = (IMG - min(IMG,[],'all'))./(max(IMG,[],'all') - min(IMG,[],'all'));

% Mostly vectorized version of above code. It is less efficient.
% % Sum full slices + part of previous slice + part of next slice
% for y = 1:size(ifull,2)
%     IMG(:,:,y) = sum(img(:,:,ifull(1,y):ifull(2,y)), 3);
% end
% % Add partial slice data (excluding index 0)
% IMG(:,:,2:end) = IMG(:,:,2:end) + img(:,:,ifull(1,2:end)-1).*reshape(ipart(1,2:end), [1,1,size(ipart(1,2:end),2)]);
% IMG = IMG + img(:,:,ifull(2,:)+1).*reshape(ipart(2,:), [1,1,size(ipart,2)]);
end
