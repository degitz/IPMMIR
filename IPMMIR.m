%% Import & Reassemble data
%%% Key
% rs: re-sized
% fs: full-sized
% og: original
%%%
%%%%%% Uses binary mask for IOU & dice, non-binary mask for RMSE
%%%%%% Transforms image & finds mask, then converts mask to 2D. Outputs only mask

clear; clc; close all;

% Parameters
tissue = 'oregano';
res = 0.069447; % resolution, mm
fov = 23; % field of view, mm
thresh = 0.25; % uCT mask threshold
roipad = 20; % roi padding
ipmethod = 'linear';  % interpolation method
rsval = [64, 128]; % resize value
a = [30, 6, 1]; % rotation resolutions & ranges for transformation
ar = [360, 60, 10];
numE = [5, 3, 1]; % number of local extrema

%%% Import uCT data
i_og_uct = niftiread("C:\Users\jacob\OneDrive - Texas A&M University\TAMU\NMR RF Lab\_Project GRMD\uCT Data\SlicerSaves\Volumes\"+tissue+"_vol.nii");
% Convert data type to double
i_og_uct = double(i_og_uct);
% Reshape & rotate array to correct format (for some reason, it outputs flipped along both axes
% This moves the "coronal" slices to the sagittal position & rotates 180 along z-axis
i_og_uct = rot90(permute(i_og_uct, [3 1 2]),2); % 'up resolution' uCT image (in reality, just normal resolution)
% Ensure data is a cube of 23mm
for i = 1:3
    arrpad = [0 0 0];
    if size(i_og_uct, i) < round(fov/res) % number of pixels for 23mm FOV of 69.447um resolution
        diff = round(fov/res) - size(i_og_uct, i);
        if rem(diff,2) == 1
            diff = diff - 1; 
            arrpad(i) = 1;
            i_og_uct = padarray(i_og_uct, arrpad, 0, 'pre');
        end
        arrpad(i) = round(diff/2);
        i_og_uct = padarray(i_og_uct, arrpad, 0, 'both');
    end
end

%%% Load MRI data
[im_name, m_name, ~, se_offset] = grmdfs(tissue, 'ssems');
[~, i_fs_mri, P, ~, ~] = loadfidfancy(im_name);
% Import mask
load(m_name); m_fs_mri = m; 
% Circular shift (so it aligns w/ GRE)
i_fs_mri = circshift(i_fs_mri, se_offset, 1);
m_fs_mri = circshift(m_fs_mri, se_offset, 1);
% Apply mask
i_fs_mri = abs(i_fs_mri.*m_fs_mri);
% Translate to center
[~, ~, diff] = CentShift(squeeze(sum(m_fs_mri,3)), squeeze(sum(i_fs_mri)));
i_fs_mri = circshift(i_fs_mri, [diff 0]); m_fs_mri = circshift(m_fs_mri, [diff 0]);

% clear arrpad fov i im_name m_name res m

%% Create binary mask for uCT
% Downscale data
% uct_im_nt = downscale3D(uct_im_fr, size(mri_im_3D)); % 'not transformed' uCT image
i_fs_uct = Norm01(imresize3(i_og_uct, [size(i_fs_mri,1) size(i_fs_mri,1) size(i_fs_mri,1)]));

% Create basic mask
m_fs_uct = i_fs_uct > thresh;

% Apply to uCT image
i_fs_uct = i_fs_uct.*m_fs_uct;

% Identify each connected region in image
CC = bwconncomp(i_fs_uct);

% Isolates largest region
[~,idx] = max(cellfun(@numel,CC.PixelIdxList)); 

% Create binary mask of only that region
m_fs_uct = false(size(i_fs_uct));
m_fs_uct(CC.PixelIdxList{idx}) = true;

% Fill holes
m_fs_uct = imfill(m_fs_uct, 'holes');

% Apply to uCT image
i_fs_uct = i_fs_uct.*m_fs_uct;

% Also center uCT data in 3D space
[m_fs_uct, i_fs_uct, ~] = CentShift(m_fs_uct, i_fs_uct);

% clear CC diff idx tr numPixels

%% Trim ROI
m_fs_mri = ROITrim(m_fs_mri, roipad);
m_fs_uct = ROITrim(m_fs_uct, roipad);

% Check if same size (likely won't be)
if size(m_fs_mri,1) > size(m_fs_uct,1)
    roidiff = size(m_fs_mri,1) - size(m_fs_uct,1);

    % Pad array
    m_fs_uct = padarray(m_fs_uct, [round(roidiff/2) round(roidiff/2) round(roidiff/2)], 0, 'both');

    % Trim any excess
    if rem(roidiff, 2) > 0
        m_fs_uct = m_fs_uct(2:end,2:end,2:end);
    end
elseif size(m_fs_mri,1) < size(m_fs_uct,1)
    roidiff = size(m_fs_uct,1) - size(m_fs_mri,1);

    % Pad array
    m_fs_mri = padarray(m_fs_mri, [round(roidiff/2) round(roidiff/2) 0], 0, 'both');

    % Trim any excess
    if rem(roidiff, 2) > 0
        m_fs_mri = m_fs_mri(2:end,2:end,:);
    end
end

% Apply trim to images
i_fs_mri = ROITrim(i_fs_mri, size(m_fs_mri));
i_fs_uct = ROITrim(i_fs_uct, size(m_fs_uct));

% clear roidiff

%% Transform & evaluate

% Resize
m_2D_mri = Norm01(sum(imresize3(m_fs_mri, [rsval(1) rsval(1) size(i_fs_mri,3)]),3));
i_rs_uct = Norm01(imresize3(i_fs_uct, [rsval(1) rsval(1) rsval(1)]));


% m_2D_mri = gpuArray(m_2D_mri); i_rs_uct = gpuArray(i_rs_uct);


% Obtain rotation angles
r = 0:a(1):(ar(1)-a(1));
Rx = repmat(r, 1, length(r)*length(r)); Rx = sort(Rx);
Ry = repmat(r, 1, length(r)); Ry = sort(Ry); Ry = repmat(Ry, 1, length(r));
Rz = repmat(r, 1, length(r)*length(r));
R = [Rx; Ry; Rz];
R = DupAngleRmv(R);

% Preallocate data
m_2D_uct = zeros(size(i_rs_uct,1), size(i_rs_uct,2), size(R,2));

% Create spatial referencing information
Rin = imref3d(size(i_rs_uct));


tic


% Transform uCT data
for i = 1:size(R,2)

    % Rotate image
    tform = affinetform3d(rigidtform3d(squeeze(R(:,i).'), [0 0 0]));
    % i_temp = imwarpcustom(i_rs_uct, Rin, tform, ipmethod, affineOutputView(size(i_rs_uct), tform, "BoundsStyle","CenterOutput"));
    i_temp = imwarp(i_rs_uct, Rin, tform, ipmethod, "OutputView", affineOutputView(size(i_rs_uct), tform, "BoundsStyle","CenterOutput"));

    % Trim to correct size (if applicable)
    if any(size(i_temp) ~= size(i_rs_uct)); i_temp = ROITrim(i_temp, size(i_rs_uct)); end

    % % Normalize, convert to 2D & add into array
    % i_temp = Norm01(i_temp);
    % i_temp = i_temp.*(i_temp > thresh);
    % m_2D_uct(:,:,i) = Norm01(sum(i_temp > thresh,3));

    % Create basic mask
    m_temp = (i_temp > thresh);

    % Fill holes
    m_temp = imfill(m_temp, 'holes');

    % Convert to 2D & add into array
    m_2D_uct(:,:,i) = sum(m_temp,3);
end

% Normalize
m_2D_uct = Norm01Slice(m_2D_uct);

% Evaluate rotatated FOV
[a_dce, a_iou, a_RMSE, B_dce, B_iou, B_RMSE, B_all] = TransEval(m_2D_uct, m_2D_mri, numE(1));


% [dce, iou, RMSE, dceB, iouB, RMSEB, allB] = gather(dce, iou, RMSE, dceB, iouB, RMSEB, allB);


% clear r ry rz Rx Ry Rz i_temp m_temp
clear m_2D_uct

% %% Transform & evaluate (4D dependent variables)
% % Obtain rotation angles
% r = 0:a(1):(360-a(1));
% [Rx,Ry,Rz] = meshgrid(r, r, r);
% R = zeros(3, length(r), length(r), length(r));
% R(1,:,:,:) = Rx; R(2,:,:,:) = Ry; R(3,:,:,:) = Rz;
% R = DupAngleRmv4D(R);
% 
% % Preallocate data
% i_2D_uct = zeros(size(i_rs_uct,1), size(i_rs_uct,2), size(R,2), size(R,3), size(R,4));
% m_2D_uct = zeros(size(i_rs_uct,1), size(i_rs_uct,2), size(R,2), size(R,3), size(R,4));
% 
% % Transform uCT data
% for i = 1:size(R,2)
%     for j = 1:size(R,3)
%         for k = 1:size(R,4)
%             if ~isnan(R(1,i,j,k))
%                 % Rotate image
%                 tform = rigidtform3d(squeeze(R(:,i,j,k)).', [0 0 0]);
%                 i_temp = imwarp(i_rs_uct, tform, ipmethod);
% 
%                 % Normalize & remove background
%                 i_temp = Norm01(i_temp);
%                 i_temp = i_temp.*(i_temp > thresh);
% 
%                 % Trim to correct size (if applicable)
%                 if any(size(i_temp) ~= size(i_rs_uct))
%                     i_temp = ROITrim(i_temp, size(i_rs_uct));
%                 end
% 
%                 % Convert to 2D & add into array
%                 m_2D_uct(:,:,i,j,k) = Norm01(sum(i_temp > thresh,3));
%                 i_2D_uct(:,:,i,j,k) = Norm01(sum(i_temp,3));
%             else
%                 m_2D_uct(:,:,i,j,k) = zeros(size(i_rs_uct,1), size(i_rs_uct,2));
%                 i_2D_uct(:,:,i,j,k) = zeros(size(i_rs_uct,1), size(i_rs_uct,2));
%             end
%         end
%     end
% end
% 
% % Evaluate rotatated FOV
% m_2D_mrirep = repmat(m_2D_mri, [1 1 size(m_2D_uct,3) size(m_2D_uct,3) size(m_2D_uct,3)]);
% area_i = squeeze(sum((m_2D_uct>0).*(m_2D_mrirep>0), [1 2]));
% area_t = (squeeze(sum(m_2D_uct, [1 2])) + squeeze(sum(m_2D_mrirep, [1 2]))) > 0;
% % IOU
% iou = area_i./(area_t - area_i);
% % dice
% dce = 2.*area_i./area_t;
% % RMSE
% RMSE = sqrt(squeeze(sum((m_2D_mrirep - m_2D_uct).^2,[1 2]))./(size(m_2D_uct,1)*size(m_2D_uct,2)));
% 
% % Find best FOV
% dceB = imregionalmax(dce, 26); % best dice values
% iouB = imregionalmax(iou, 26); % best iou values
% RMSEB = imregionalmin(RMSE, 26); % best RMSE values
% allB = dceB.*iouB.*RMSEB; allB = find(allB == 1);
% 
% clear r ry rz Rx Ry Rz i_temp m_temp

%% Secondary Transformations
% medium size image analysis, viewing high resolution rotations around local minima

% Second resizing
m_2D_mri = Norm01(sum(imresize3(m_fs_mri, [rsval(2) rsval(2) size(i_fs_mri,3)]),3));
i_rs_uct = Norm01(imresize3(i_fs_uct, [rsval(2) rsval(2) rsval(2)]));

for m = 1:length(B_all)
    % Obtain rotation angles
    rx = ((R(1,B_all(m))-ar(2)/2):a(2):(R(1,B_all(m))+ar(2)/2)); ry = ((R(2,B_all(m))-ar(2)/2):a(2):(R(2,B_all(m))+ar(2)/2)); rz = ((R(3,B_all(m))-ar(2)/2):a(2):(R(3,B_all(m))+ar(2)/2));
    Rx = repmat(rx, 1, length(ry)*length(rz)); Rx = sort(Rx);
    Ry = repmat(ry, 1, length(rz)); Ry = sort(Ry); Ry = repmat(Ry, 1, length(ry));
    Rz = repmat(rz, 1, length(rx)*length(ry));

    % Add to old angles & remove duplicates
    Rold = R;
    R = DupAngleRmv(horzcat(Rold, [Rx; Ry; Rz]));

    % Isolate new angles
    Rnew = R(:, (size(Rold,2)+1):end);
    Rnew = horzcat(R(:,B_all(m)), Rnew);

    % Preallocate data
    m_2D_uct_temp = zeros(size(i_rs_uct,1), size(i_rs_uct,2), size(Rnew,2), 'like', i_rs_uct);

    % Create spatial referencing information
    Rin = imref3d(size(i_rs_uct));

    % Transform uCT data
    for i = 1:size(Rnew,2)

        % Rotate image
        tform = affinetform3d(rigidtform3d(squeeze(Rnew(:,i).'), [0 0 0]));
        % i_temp = imwarpcustom(i_rs_uct, Rin, tform, ipmethod, affineOutputView(size(i_rs_uct), tform, "BoundsStyle","CenterOutput"));
        i_temp = imwarp(i_rs_uct, Rin, tform, ipmethod, "OutputView", affineOutputView(size(i_rs_uct), tform, "BoundsStyle","CenterOutput"));

        % Trim to correct size (if applicable)
        if any(size(i_temp) ~= size(i_rs_uct)); i_temp = ROITrim(i_temp, size(i_rs_uct)); end
        
        % % Normalize, convert to 2D & add into array
        % i_temp = Norm01(i_temp);
        % i_temp = i_temp.*(i_temp > thresh);
        % m_2D_uct(:,:,i) = Norm01(sum(i_temp > thresh,3));

        % Create basic mask
        m_temp = (i_temp > thresh);

        % Fill holes
        m_temp = imfill(m_temp, 'holes');

        % Convert to 2D & add into array
        m_2D_uct_temp(:,:,i) = sum(m_temp,3);
    end

    if m == 1
        m_2D_uct = m_2D_uct_temp;
    else
        % Add images to array
        m_2D_uct = cat(3, m_2D_uct, m_2D_uct_temp(:,:,2:end));
    end
end

% Normalize
m_2D_uct = Norm01Slice(m_2D_uct);

% Evaluate rotatated FOV
[a_dce, a_iou, a_RMSE, B_dce, B_iou, B_RMSE, B_all] = TransEval(m_2D_uct, m_2D_mri, numE(2), a_dce, a_iou, a_RMSE);

% clear rx ry rz Rx Ry Rz RB i_2D_uct_temp m_2D_uct_temp i_temp m_temp Rold Rnew

%% Tertiary Transformations
% Full size image analysis, viewing high resolution rotations around local minima

for m = 1:length(B_all)
    % Obtain rotation angles
    rx = ((R(1,B_all(m))-ar(3)/2):a(3):(R(1,B_all(m))+ar(3)/2)); ry = ((R(2,B_all(m))-ar(3)/2):a(3):(R(2,B_all(m))+ar(3)/2)); rz = ((R(3,B_all(m))-ar(3)/2):a(3):(R(3,B_all(m))+ar(3)/2));
    Rx = repmat(rx, 1, length(ry)*length(rz)); Rx = sort(Rx);
    Ry = repmat(ry, 1, length(rz)); Ry = sort(Ry); Ry = repmat(Ry, 1, length(ry));
    Rz = repmat(rz, 1, length(rx)*length(ry));

    % Add to old angles & remove duplicates
    Rold = R;
    R = DupAngleRmv(horzcat(Rold, [Rx; Ry; Rz]));

    % Isolate new angles
    Rnew = R(:, (size(Rold,2)+1):end);
    Rnew = horzcat(R(:,B_all(m)), Rnew);

    % Preallocate data
    m_2D_uct_temp = zeros(size(i_fs_uct,1), size(i_fs_uct,2), size(Rnew,2), 'like', i_fs_uct);

    % Create spatial referencing information
    Rin = imref3d(size(i_fs_uct));

    % Transform uCT data
    for i = 1:size(Rnew,2)

        % Rotate image
        tform = affinetform3d(rigidtform3d(squeeze(Rnew(:,i).'), [0 0 0]));
        % i_temp = imwarpcustom(i_fs_uct, Rin, tform, ipmethod, affineOutputView(size(i_fs_uct), tform, "BoundsStyle","CenterOutput"));
        i_temp = imwarp(i_fs_uct, Rin, tform, ipmethod, "OutputView", affineOutputView(size(i_fs_uct), tform, "BoundsStyle","CenterOutput"));

        % Trim to correct size (if applicable)
        if any(size(i_temp) ~= size(i_fs_uct)); i_temp = ROITrim(i_temp, size(i_fs_uct)); end

        % % Normalize, convert to 2D & add into array
        % i_temp = Norm01(i_temp);
        % i_temp = i_temp.*(i_temp > thresh);
        % m_2D_uct(:,:,i) = Norm01(sum(i_temp > thresh,3));

        % Create basic mask
        m_temp = (i_temp > thresh);

        % Fill holes
        m_temp = imfill(m_temp, 'holes');

        % Convert to 2D & add into array
        m_2D_uct_temp(:,:,i) = sum(m_temp,3);
    end

    if m == 1
        m_2D_uct = m_2D_uct_temp;
    else
        % Add images to array
        m_2D_uct = cat(3, m_2D_uct, m_2D_uct_temp(:,:,2:end));
    end
end


toc


% Normalize
m_2D_uct = Norm01Slice(m_2D_uct);
m_2D_mri = Norm01(sum(m_fs_mri,3));

% Evaluate rotatated FOV
[a_dce, a_iou, a_RMSE, B_dce, B_iou, B_RMSE, B] = TransEval(m_2D_uct, m_2D_mri, numE(3), a_dce, a_iou, a_RMSE);

% clear rx ry rz Rx Ry Rz RB i_2D_uct_temp m_2D_uct_temp i_temp m_temp Rold Rnew

%% Show best FOV

data = table([1:length(a_iou)].', R.', a_dce, a_iou, a_RMSE); %#ok<NBRAK1>
data.Properties.VariableNames = {'Num', 'Rotation', 'dice', 'iou', 'RMSE'};

% take data off gpu for non-linear interpolation
if ~strcmp(ipmethod,'linear') && isgpuarray(i_fs_uct)
    i_fs_uct = gather(i_fs_uct);
end

% Create MRI & uCT images for display
tform = affinetform3d(rigidtform3d(squeeze(R(:,B).'), [0 0 0]));
i_final_uct = imwarpcustom(i_fs_uct, imref3d(size(i_fs_uct)), tform, 'cubic', affineOutputView(size(i_fs_uct), tform, "BoundsStyle","CenterOutput"));
i_final_uct = Norm01(i_final_uct); i_final_uct = i_final_uct.*(i_final_uct > thresh);
if any(size(i_final_uct) ~= size(i_fs_uct)); i_final_uct = ROITrim(i_final_uct, size(i_fs_uct)); end
i_ds_uct = Norm01(sum(i_final_uct,3));
i_ds_mri = Norm01(sum(i_fs_mri, 3));

% Show figures
figure('Name', 'MRI_mb', 'NumberTitle', 'off'); imagesc(m_2D_mri>0); axis square; axis off; colormap gray; title('2D MRI Mask (binary)');
figure('Name', 'uCT_mb', 'NumberTitle', 'off'); imagesc(m_2D_uct(:,:,B)>0); axis square; axis off; colormap gray; title('2D uCT Mask (binary)');
figure('Name', 'com_mb', 'NumberTitle', 'off'); imagesc(imfuse(m_2D_mri>0, m_2D_uct(:,:,B)>0)); axis square; axis off; title('Image Overlap (binary)'); hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'MRI Image','uCT Image'}); hold off;
figure('Name', 'MRI_m', 'NumberTitle', 'off'); imagesc(m_2D_mri); axis square; axis off; colormap gray; title('2D MRI Mask');
figure('Name', 'uCT_m', 'NumberTitle', 'off'); imagesc(m_2D_uct(:,:,B)); axis square; axis off; colormap gray; title('2D uCT Mask');
figure('Name', 'com_m', 'NumberTitle', 'off'); imagesc(imfuse(m_2D_mri, m_2D_uct(:,:,B))); axis square; axis off; title('Image Overlap'); hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'MRI Image','uCT Image'}); hold off;
figure('Name', 'MRI_i', 'NumberTitle', 'off'); imagesc(i_ds_mri); axis square; axis off; colormap gray; title('2D MRI');
figure('Name', 'uCT_i', 'NumberTitle', 'off'); imagesc(i_ds_uct); axis square; axis off; colormap gray; title('2D uCT');
figure('Name', 'com_i', 'NumberTitle', 'off'); imagesc(imfuse(i_ds_mri, i_ds_uct)); axis square; axis off; title('Image Overlap'); hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'MRI Image','uCT Image'}); hold off;

% dce_scale = a_dce.^4; dce_scale = (dce_scale./max(dce_scale)).*40;
% iou_scale = a_iou.^4; iou_scale = (iou_scale./max(iou_scale)).*40;
% RMSE_scale = a_RMSE.^4; RMSE_scale = (RMSE_scale./min(RMSE_scale)).*40;
% figure('Name', 'dice', 'NumberTitle', 'off'); scatter3(Rx,Ry,Rz,dce_scale, a_dce, 'filled'); xlim([0 360]); ylim([0 360]); zlim([0 360]); xlabel('X Rotation (Degrees)'); ylabel('Y Rotation (Degrees)'); zlabel('Z Rotation (Degrees)'); cb = colorbar; cb.Label.String = 'dice'; colormap(flipud(turbo));
% figure('Name', 'iou', 'NumberTitle', 'off'); scatter3(Rx,Ry,Rz,iou_scale, a_iou, 'filled'); xlim([0 360]); ylim([0 360]); zlim([0 360]); xlabel('X Rotation (Degrees)'); ylabel('Y Rotation (Degrees)'); zlabel('Z Rotation (Degrees)'); cb = colorbar; cb.Label.String = 'iou'; colormap(flipud(turbo));
% figure('Name', 'RMSE', 'NumberTitle', 'off'); scatter3(Rx,Ry,Rz,RMSE_scale, a_RMSE, 'filled'); xlim([0 360]); ylim([0 360]); zlim([0 360]); xlabel('X Rotation (Degrees)'); ylabel('Y Rotation (Degrees)'); zlabel('Z Rotation (Degrees)'); cb = colorbar; cb.Label.String = 'RMSE'; colormap(turbo);

% clear i j k i2 j2 m diff tissue CC idx numPixels i_temp m_temp tform im_name m_name m

%% Gradient
% [FX, FY, FZ] = gradient(iou, [a(1) a(1) a(1)]);
%% Functions
% Assuming both inputs are same FOV, decreases resolution of first image to match that of second image
function img_newres = Downscale3D(img_oldres, imsize)
dimlen = size(img_oldres); % find image size

% Check if input is square
if min(dimlen) ~= max(dimlen)

    % Find size difference
    diff = max(dimlen) - min(dimlen); diff1 = floor(diff/2);

    % Check if difference is odd (if so, remove extra from one side)
    extra = rem(diff, 2);

    % Trim FOV
    [~, Mdim] = max(dimlen);
    if Mdim == 1
        if dimlen(2) == dimlen(1)
            img_oldres = img_oldres(1+diff1:dimlen-diff1-extra, 1+diff1:dimlen-diff1-extra, :);
        elseif dimlen(3) == dimlen(1)
            img_oldres = img_oldres(1+diff1:dimlen-diff1-extra, :, 1+diff1:dimlen-diff1-extra);
        else
            img_oldres = img_oldres(1+diff1:dimlen-diff1-extra, :, :);
        end
    elseif Mdim == 2
        if dimlen(1) == dimlen(2)
            img_oldres = img_oldres(1+diff1:dimlen-diff1-extra, 1+diff1:dimlen-diff1-extra, :);
        elseif dimlen(3) == dimlen(2)
            img_oldres = img_oldres(:, 1+diff1:dimlen-diff1-extra, 1+diff1:dimlen-diff1-extra);
        else
            img_oldres = img_oldres(:, 1+diff1:dimlen-diff1-extra, :);
        end
    else % Mdim ==
        if dimlen(1) == dimlen(3)
            img_oldres = img_oldres(1+diff1:dimlen-diff1-extra, :, 1+diff1:dimlen-diff1-extra);
        elseif dimlen(2) == dimlen(3)
            img_oldres = img_oldres(:, 1+diff1:dimlen-diff1-extra, 1+diff1:dimlen-diff1-extra);
        else
            img_oldres = img_oldres(:, :, 1+diff1:dimlen-diff1-extra);
        end
    end
end

dimlen = size(img_oldres,1); % find image size

kspace = ifftshift(fftn(ifftshift(img_oldres))); % kspace of og data

% Find size difference
diff = size(img_oldres, 1) - imsize(1); diff1 = floor(diff/2);

% Check if difference is odd (if so, remove extra from one side)
extra = rem(diff, 2);

% Decrease resolution
kspace = kspace(1+diff1:dimlen-diff1-extra, 1+diff1:dimlen-diff1-extra, 1+diff1:dimlen-diff1-extra);

% Convert from FD
img_newres = abs(fftshift(ifftn(fftshift(kspace))));
end

% Assuming both inputs are same FOV, decreases resolution of first image to match that of second image
function img_newres = Downscale2D(img_oldres, imsize)
dimlen = size(img_oldres,1); % find image size

kspace = ifftshift(fftn(ifftshift(img_oldres))); % kspace of og data

% Find size difference
diff = size(img_oldres, 1) - imsize(1); diff1 = floor(diff/2);

% Check if difference is odd (if so, remove extra from one side)
extra = rem(diff, 2);

% Decrease resolution
kspace = kspace(1+diff1:dimlen-diff1-extra, 1+diff1:dimlen-diff1-extra, :);

% Convert from FD
img_newres = abs(fftshift(ifftn(fftshift(kspace))));
end

% Centers 2D data
function [mcs, imcs, diff] = CentShift(m, im)
% Find center of ROI
if size(m, 3) > 1
    cent = regionprops3(m, 'Centroid', 'Volume'); 

    % Check if multiple ROIs are being detected
    if size(cent,1) > 1
        trueROI = find(cent.Volume == max(cent.Volume));
    else
        trueROI = 1;
    end

    cent = round(cent.Centroid(trueROI,:));

    % Find offset from center
    diff = [round(size(m,1)/2)-cent(2), round(size(m,2)/2)-cent(1), round(size(m,3)/2)-cent(3)];
else
    
    cent = regionprops(m, 'Centroid'); 
    cent = round(cat(1, cent.Centroid));

    % Find offset from center
    diff = [round(size(m,1)/2)-cent(2), round(size(m,2)/2)-cent(1)];
end

% Shift mask & image to center
mcs = circshift(m, diff); imcs = circshift(im, diff);
end

% Normalize input data to range of 0-1
function X = Norm01(x)
X = (x - min(x,[],'all'))./max(x,[],'all');
end

% Normalize input data to range of 0-1 by slice
function X = Norm01Slice(x)
X = (x - min(x,[],3))./max(x,[],3);
end

% Trim image down to specified size. Buff is either equal padding on all axis, or dimension size that roi should be padded to
function T = ROITrim(ut, buff)
% Reset ROI range
roi = [1, size(ut,1); 1, size(ut,2); 1, size(ut,3)];

% Find mask limits
[~, roi(2,1)] = find(squeeze(sum(ut, 3)), 1, 'first');
[~, roi(2,2)] = find(squeeze(sum(ut, 3)), 1, 'last');
[~, roi(1,1)] = find(squeeze(sum(ut, 3)).', 1, 'first');
[~, roi(1,2)] = find(squeeze(sum(ut, 3)).', 1, 'last');
roi(3,1) = find(squeeze(sum(ut, [1 2])), 1, 'first');
roi(3,2) = find(squeeze(sum(ut, [1 2])), 1, 'last');

% Find axis sizes
roisz = [(roi(1,2) - roi(1,1) + 1) (roi(2,2) - roi(2,1) + 1) (roi(3,2) - roi(3,1) + 1)];
[~, dmax] = max(roisz);

if max(roisz) ~= min(roisz)
    % Check dimensionality
    if size(ut,3)*5 < size(ut,1)
        % Determine which axis is smaller
        if dmax == 1
            dmin = 2;
        else % dmax == 2
            dmin = 1;
        end

        % Find difference
        roidiff = roisz(dmax) - roisz(dmin);

        % Add difference to smaller dimension
        roi(dmin,1) = roi(dmin,1) - round(roidiff/2);
        roi(dmin,2) = roi(dmin,2) + round((roidiff-1)/2);
    else
        dmin = [0 0];
        % Determine which axis(es) are smaller
        if dmax == 1
            if roisz(2) ~= max(roisz)
                dmin(1) = 2;
                if roisz(3) ~= max(roisz)
                    dmin(2) = 3;
                else
                    dmin(2) = -3;
                end
            elseif roisz(3) ~= max(roisz)
                dmin(1) = 3; dmin(2) = -2;
            end
        elseif dmax == 2
            dmin = [1 -3];
            if roisz(3) ~= max(roisz)
                dmin(2) = 3;
            end
        else % dmax == 3
            dmin = [1 2];
        end

        % Find first difference
        roidiff = roisz(dmax) - roisz(dmin(1));

        % Add difference to smaller dimension
        roi(dmin(1),1) = roi(dmin(1),1) - round(roidiff/2);
        roi(dmin(1),2) = roi(dmin(1),2) + round((roidiff-1)/2);

        % Update axis sizes
        roisz = [(roi(1,2) - roi(1,1) + 1) (roi(2,2) - roi(2,1) + 1) (roi(3,2) - roi(3,1) + 1)];

        % Check if only one smaller axis
        if dmin(2) > 0
            % Find second difference
            roidiff = roisz(dmax) - roisz(dmin(2));

            % Add difference to smaller dimension
            roi(dmin(2),1) = roi(dmin(2),1) - round(roidiff/2);
            roi(dmin(2),2) = roi(dmin(2),2) + round((roidiff-1)/2);
        end
    end
end

% Update axis sizes
roisz = [(roi(1,2) - roi(1,1) + 1) (roi(2,2) - roi(2,1) + 1) (roi(3,2) - roi(3,1) + 1)];

if size(ut,3)*5 < size(ut,1)
    % Pad ROI so tissue is not on edge of figures
    [t1, t2] = size(buff);
    if (t1 == 1) && (t2 == 1)
        roi(1:2,1) = roi(1:2,1) - buff;
        roi(1:2,2) = roi(1:2,2) + buff;

    elseif any((buff - roisz) < 0)
        error('Rotated ROI is larger than expected array size. Please increase padding on initial ROI trim.')
    else
        roidiff = buff - roisz;

        % Change to correct size (& take into account odd difference)
        roi(1:2,1) = roi(1:2,1) - round((roidiff(1:2)-1)./2).';
        roi(1:2,2) = roi(1:2,2) + round(roidiff(1:2)./2).';
    end
else
    % Pad ROI so tissue is not on edge of figures
    [t1, t2] = size(buff);
    if (t1 == 1) && (t2 == 1)
        roi(:,1) = roi(:,1) - buff;
        roi(:,2) = roi(:,2) + buff;

    elseif any((buff - roisz) < 0)
        error('ROI is larger than expected array size. Please increase padding on initial ROI trim.')
    else
        roidiff = buff - roisz;

        % Change to correct size (& take into account odd difference)
        roi(:,1) = roi(:,1) - round((roidiff-1)./2).';
        roi(:,2) = roi(:,2) + round(roidiff./2).';
    end
end

% Check for values out of bounds
if any(roi(:,1) <= 0)
    x = find(roi(:,1) <= 0);
    roi(x,:) = roi(x,:) + 1;
elseif any(roi(:,2) > size(ut).')
    x = find(roi(:,2) > size(ut).');
    roi(x,:) = roi(x,:) - 1;
end

% Trim data
T = ut(roi(1,1):roi(1,2), roi(2,1):roi(2,2), roi(3,1):roi(3,2));
end

% Remove angles that would produce duplicate rotational matrices
function Rcorr = DupAngleRmv(R)
Rcorr = zeros(size(R), 'like', R); I = 1;
eulcheck = zeros(size(R,1), size(R,1), size(R,2), 'like', R);
RM = zeros(size(R,1), size(R,1), size(R,2), 'like', R);

% Calculate rotational matrices
RM(1,:,:) = [cosd(R(2,:)).*cosd(R(3,:)); sind(R(1,:)).*sind(R(2,:)).*cosd(R(3,:))-cosd(R(1,:)).*sind(R(3,:)); cosd(R(1,:)).*sind(R(2,:)).*cosd(R(3,:))+sind(R(1,:)).*sind(R(3,:))];
RM(2,:,:) = [cosd(R(2,:)).*sind(R(3,:)); sind(R(1,:)).*sind(R(2,:)).*sind(R(3,:))+cosd(R(1,:)).*cosd(R(3,:)); cosd(R(1,:)).*sind(R(2,:)).*sind(R(3,:))-sind(R(1,:)).*cosd(R(3,:))];
RM(3,:,:) = [-sind(R(2,:)); sind(R(1,:)).*cosd(R(2,:)); cosd(R(1,:)).*cosd(R(2,:))];

% Correct for any rounding errors between GPU & CPU
RM = round(RM*10e9);

% Check for duplicates
for i = 2:size(R,2)
    % Check if rotational matrix is equal to any others already included
    if all(sum(abs(eulcheck - RM(:,:,i)), [1 2]) ~= 0, 'all') && all(sum(abs(eulcheck + RM(:,:,i)), [1 2]) ~= 0, 'all') 
        I = I + 1;
        eulcheck(:,:,i) = RM(:,:,i);
        Rcorr(:,I) = R(:,i);
    end
end
Rcorr = Rcorr(:,1:I);
end

% Remove angles that would produce duplicate rotational matrices (4D data)
function Rcorr = DupAngleRmv4D(R)
% Check for unique rotational matrices
Rcorr = zeros(size(R), 'like', R);
eulcheck = zeros(size(R,1), size(R,1), size(R,2), size(R,3), size(R,4), 'like', R);
RM = zeros(size(R,1), size(R,1), size(R,2), size(R,3), size(R,4), 'like', R);

RM(1,:,:,:,:) = [cosd(R(2,:,:,:)).*cosd(R(3,:,:,:)); sind(R(1,:,:,:)).*sind(R(2,:,:,:)).*cosd(R(3,:,:,:))-cosd(R(1,:,:,:)).*sind(R(3,:,:,:)); cosd(R(1,:,:,:)).*sind(R(2,:,:,:)).*cosd(R(3,:,:,:))+sind(R(1,:,:,:)).*sind(R(3,:,:,:))];
RM(2,:,:,:,:) = [cosd(R(2,:,:,:)).*sind(R(3,:,:,:)); sind(R(1,:,:,:)).*sind(R(2,:,:,:)).*sind(R(3,:,:,:))+cosd(R(1,:,:,:)).*cosd(R(3,:,:,:)); cosd(R(1,:,:,:)).*sind(R(2,:,:,:)).*sind(R(3,:,:,:))-sind(R(1,:,:,:)).*cosd(R(3,:,:,:))];
RM(3,:,:,:,:) = [-sind(R(2,:,:,:)); sind(R(1,:,:,:)).*cosd(R(2,:,:,:)); cosd(R(1,:,:,:)).*cosd(R(2,:,:,:))];

for i = 2:size(R,2)
    for j = 1:size(R,3)
        for k = 1:size(R,4)
            % Check if rotational matrix is equal to any others already included
            if all(sum(abs(eulcheck - RM(:,:,i,j,k)), [1 2]) ~= 0, 'all') && all(sum(abs(eulcheck + RM(:,:,i,j,k)), [1 2]) ~= 0, 'all')
                eulcheck(:,:,i,j,k) = RM(:,:,i,j,k);
                Rcorr(:,i,j,k) = R(:,i,j,k);
            else
                Rcorr(:,i,j,k) = NaN(3,1);
            end
        end
    end
end
end

% Evaluate transformations
function [a_dce, a_iou, a_RMSE, B_dce, B_iou, B_RMSE, B_all] = TransEval(i_tran, i_true, numE, a_dce, a_iou, a_RMSE)
% Check if first use
if nargin == 3
    a_dce = []; a_iou = []; a_RMSE = []; i = 1;
else
    i = 2;
end
    
% Evaluate rotatated FOV
i_true = repmat(i_true, [1 1 size(i_tran,3)]);
area_i = squeeze(sum((i_tran>0).*(i_true>0), [1 2]));
area_t = (squeeze(sum(i_tran>0, [1 2])) + squeeze(sum(i_true>0, [1 2])));
% IOU
ioutemp = area_i./(area_t - area_i); a_iou = [a_iou; ioutemp(i:end)];
% dice
dcetemp = 2.*area_i./area_t; a_dce = [a_dce; dcetemp(i:end)];
% RMSE
RMSEtemp = sqrt(squeeze(sum((i_true - i_tran).^2,[1 2]))./(size(i_tran,1)*size(i_tran,2))); a_RMSE = [a_RMSE; RMSEtemp(i:end)];

% Find best FOV
[~, B_dce] = maxk(a_dce, 60); % best dice values
[~, B_iou] = maxk(a_iou, 60); % best iou values
[~, B_RMSE] = mink(a_RMSE, 60); % best RMSE values
B_all = intersect(B_dce, B_iou); B_all = intersect(B_all, B_RMSE);

% Check for no intersections
if isempty(B_all)
    B_RMSE = B_RMSE(1:30);
    [~, B_all] = maxk(a_iou(B_RMSE), numE);
elseif length(B_all) > numE
    B_all = B_all(1:numE); 
end
end

% Custom version of 'imwarp()'; more efficient & works better w/ GPU
function outputImage = imwarpcustom(inputImage, Rin, tform, method, Rout)

% Define affine transformation that maps from intrinsic system of output image to world system of output image.
tIntrinsictoWorldOutput = Rout.TransformIntrinsicToWorld;

% Define affine transformation that maps from world system of input image to intrinsic system of input image.
tWorldToIntrinsicInput = Rin.TransformWorldToIntrinsic;

% Form transformation to go from output intrinsic to input intrinsic.
% NOTE: tComp is in post-multiply form.
tComp = tIntrinsictoWorldOutput/tform.A'*tWorldToIntrinsicInput; %tComp = Rout.TransformIntrinsicToWorld/tform.A'*Rin.TransformWorldToIntrinsic;
tComp(:,4)=[0 0 0 1]; % avoid round off issues due to inversion above

% Find the transform that takes from input intrinsic to output intrinsic
tformComposite = affinetform3d(tComp');

% Form plaid grid of intrinsic points in output image.
[dstXIntrinsic,dstYIntrinsic,dstZIntrinsic] = meshgrid(1:Rout.ImageSize(2), 1:Rout.ImageSize(1), 1:Rout.ImageSize(3));
[srcXIntrinsic,srcYIntrinsic, srcZIntrinsic] = transformPointsForward(tformComposite, dstXIntrinsic,dstYIntrinsic, dstZIntrinsic);

% Interpolate to produce transformed image
outputImage = interp3(inputImage,srcXIntrinsic,srcYIntrinsic,srcZIntrinsic,method);

% Replace NaNs with zeros
outputImage(isnan(outputImage)) = 0;

% outputImage = images.internal.interp3d(inputImage,srcXIntrinsic,srcYIntrinsic,srcZIntrinsic,method,fillValues, SmoothEdges);
% works in a similar manner as interp3 when smoothing is turned off except it automatically replaces NaNs w/ zeros
end
