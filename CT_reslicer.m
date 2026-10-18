function [i_CT, m_CT, Bpix, Bofs] = CT_reslicer(i_ct, m_ct, i_mri, m_mri, r, etype)
% ct - ct data
% mri - mri data
% r - number of ct slices per mri slice (mri res / ct res)

% Select evaluation type
% Select evaluation metric
switch etype
    case 'Dice'
        SliceEval = @(v,w,x,y,z) dicecustom(v,w,x,y,z);
    case 'iRMSE'
        SliceEval = @(v,w,x,y,z) irmsecustom(v,w,x,y,z);
    case 'IoU'
        SliceEval = @(v,w,x,y,z) ioucustom(v,w,x,y,z);
    case'Dice_iRMSE'
        SliceEval = @(v,w,x,y,z) dicecustom(v,w,x,y,z)+irmsecustom(v,w,x,y,z);
end

% Pad ct data to give us room to play
i_ct = padarray(i_ct, [0 0 ceil(r)], 0, 'both');
m_ct = padarray(m_ct, [0 0 ceil(r)], false, 'both');

% Determine indices of start/stop for mri slices in ct
ibase = 0:r:size(i_ct,3); % base starts & stops
ifull(1,:) = ceil(ibase(1:end-1))+1; % first full slices
ifull(2,:) = floor(ibase(2:end)); % last full slices

% Find non-integer amounts contributing to slices
ipart(1,:) = ceil(ibase(1:end-1))-ibase(1:end-1);
ipart(2,:) = ibase(2:end)-floor(ibase(2:end));

% Preallocate data
met = zeros(size(i_ct,3)-ifull(2,end) - 1, 1);
Bofs = zeros(size(i_ct,3)-ifull(2,end), 1);

% Compare CT slice assignments with available subslice offsets. Last offset
% is excluded, since partial pixel addition would extend past matrix bounds.
for pix = 0:size(i_ct,3)-ifull(2,end)-1
    % fprintf('\nTesting a subslice offset of %i pixel(s).\n', pix)
    i_CT = SliceCombine(i_ct, ifull+pix, ipart);
    m_CT = SliceCombine(m_ct, ifull+pix, ipart);

    % For extra slices in CT data, determine offset with MRI
    % Preallocate data
    met_ofs = zeros(size(i_CT,3)-size(i_mri,3)+1,1);
    for slc = 0:size(i_CT,3)-size(i_mri,3)

        % Isolate current number of slices
        i_CT_ofs = i_CT(:,:,1+slc:end-(size(i_CT,3)-size(i_mri,3))+slc);
        m_CT_ofs = m_CT(:,:,1+slc:end-(size(m_CT,3)-size(i_mri,3))+slc);

        met_ofs(slc+1) = SliceEval(i_CT_ofs, m_CT_ofs, i_mri, m_mri, 3);
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

    % Determine best
    Bofs(pix+1) = find(met_ofs == max(met_ofs), 1) - 1;
    i_CT = i_CT(:,:,1+Bofs(pix+1):end-(size(i_CT,3)-size(i_mri,3))+Bofs(pix+1));
    m_CT = m_CT(:,:,1+Bofs(pix+1):end-(size(m_CT,3)-size(i_mri,3))+Bofs(pix+1));

    % fprintf('For a %i pixel subslice offset, a %i slice offset was best.\n', pix, Bofs(pix+1))

    % Evaluate
    met(pix+1) = SliceEval(i_CT, m_CT, i_mri, m_mri, 3);
end

% Isolate best
Bpix = find(met == max(met));
Bofs = Bofs(Bpix); Bpix = Bpix-1;

% Isolate correct slice & pixel values
firstslice = 1 + Bofs;
lastslice = size(i_mri,3) + Bofs;
% firstpixel = ifull(1,firstslice) + Bpix;
% lastpixel = ifull(2,lastslice) + Bpix;
% fprintf('\nOverall best results obtained with offsets of %i pixel(s) and %i slice(s).\n', Bpix, Bofs)
% fprintf('This includes slices %i to %i and pixels %i to %i, plus %.4f and %.4f of pixels %i and %i, respectively.\n', firstslice, lastslice, firstpixel, lastpixel, ipart(1,firstslice), ipart(2,lastslice), firstpixel-1, lastpixel+1)
i_CT = SliceCombine(i_ct, ifull+Bpix, ipart);
m_CT = SliceCombine(m_ct, ifull+Bpix, ipart);
i_CT = i_CT(:,:,firstslice:lastslice);
m_CT = m_CT(:,:,firstslice:lastslice);
end

% Combine microCT slices into MRI slices
function IMG = SliceCombine(img, ifull, ipart)
% Preallocate data
IMG = zeros(size(img,1), size(img,2), size(ifull,2), 'like', img);

% Assign slices to CT data
IMG(:,:,1) = sum(img(:,:,ifull(1,1):ifull(2,1)), 3) + img(:,:,ifull(2,1)+1).*ipart(2,1); % hard code first iteration to avoid if statement (first slice should start at an integer)
for y = 2:size(ifull,2)
    % Sum full slices + part of previous slice + part of next slice
    IMG(:,:,y) = sum(img(:,:,ifull(1,y):ifull(2,y)), 3) + img(:,:,ifull(1,y)-1).*ipart(1,y) + img(:,:,ifull(2,y)+1).*ipart(2,y);
end

% Normalize
IMG = (IMG - min(IMG,[],'all'))./max(IMG,[],'all');

% Mostly vectorized version of above code. It is less efficient.
% % Sum full slices + part of previous slice + part of next slice
% for y = 1:size(ifull,2)
%     IMG(:,:,y) = sum(img(:,:,ifull(1,y):ifull(2,y)), 3);
% end
% % Add partial slice data (excluding index 0)
% IMG(:,:,2:end) = IMG(:,:,2:end) + img(:,:,ifull(1,2:end)-1).*reshape(ipart(1,2:end), [1,1,size(ipart(1,2:end),2)]);
% IMG = IMG + img(:,:,ifull(2,:)+1).*reshape(ipart(2,:), [1,1,size(ipart,2)]);
end
