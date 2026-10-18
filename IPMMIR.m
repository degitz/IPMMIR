%% Import & Reassemble data
%%% Key
% rs: re-sized
% fs: full-sized
% og: original
%%%
%%%%%% Uses binary mask for IOU & dice, non-binary mask for RMSE
%%%%%% Transforms image & finds mask, then converts mask to 2D. Outputs only mask
clear; clc; close all;

tic

% Parameters
tissue = 'lyle';
res = 0.069447; % resolution, mm
thresh = 0.25; % uCT mask threshold
roipad = 20; % roi padding
ipmethod = 'linear';  % interpolation method
rsval = [64, 96, 128]; % resize value
a = [30, 10, 2, 1]; % rotation resolutions & ranges for transformation
ar = [360, 60, 12, 6];
numE = [5, 4, 3, 1]; % number of local extrema
nume = 10; % number of local extrema before intersection
S = length(a); % number of iterations

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

fov = P.FOV.fov.ro*10; % field of view, mm

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

clear diff i fov m m_name im_name res arrpad P

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

% Apply to uCT image
i_fs_uct = i_fs_uct.*m_fs_uct;

clear CC idx

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

% Lastly, center data in 3D space
% Note: it is important to center the data after the ROI is trimmed. If
% done so beforehand, the ROI trimming will off-center the volume to ensure
% equal padding on all sides
[m_fs_mri, i_fs_mri, ~] = CentShift(m_fs_mri, i_fs_mri);
[m_fs_uct, i_fs_uct, ~] = CentShift(m_fs_uct, i_fs_uct);

fprintf('\nFinished loading & pre-processing data.\n')

toc

clear roidiff

%% Transform & Evaluate
Bi_all = 1; R = [180; 180; 180]; a_iou = []; a_dce = []; a_RMSE = []; fr = []; % initial

for s = 1:S
    fprintf('\n--------------------------- Phase %i ---------------------------\n', s)

    % Resize (if applicable)
    if s < S
        m_test_mri = Norm01(sum(imresize3(m_fs_mri, [rsval(s) rsval(s) size(m_fs_mri,3)]),3));
        i_test_uct = Norm01(imresize3(i_fs_uct, [rsval(s) rsval(s) rsval(s)]));
    else
        m_test_mri = Norm01(sum(m_fs_mri,3));
        i_test_uct = i_fs_uct;
    end

    % Create spatial referencing information
    Rin = imref3d(size(i_test_uct));

    fprintf('\nTesting the following degree range(s) with a %i degree step size:\n', a(s))

    % Determine if any groups overlap
    ccheck = false(length(Bi_all), 1);
    Rbase = R(:,Bi_all); % angle base
    Rmax = Rbase + ar(s)/2; % angle max
    Rmin = Rbase - ar(s)/2; % angle min

    % Correct for 360 values
    Rmax(Rmax==360) = 0;
    r = zeros(6, length(Bi_all));

    % Find all rotation angles with no repeats
    for m = 1:length(Bi_all)
        % Obtain rotation angles
        rx = (R(1,Bi_all(m))-ar(s)/2):a(s):(R(1,Bi_all(m))+ar(s)/2);
        ry = (R(2,Bi_all(m))-ar(s)/2):a(s):(R(2,Bi_all(m))+ar(s)/2);
        rz = (R(3,Bi_all(m))-ar(s)/2):a(s):(R(3,Bi_all(m))+ar(s)/2);

        % Skips for overlapping groups
        if ccheck(m) == true
            continue
        end

        if m < length(Bi_all)
            for l = 1:(length(Bi_all)-m)
                % Skip checks with already combined groups
                if ccheck(m+l) == true
                    continue
                else
                    % Create testing variables
                    % Find equal angles
                    Tequ = Rbase(:,m+l) == Rbase(:,m);
                    % Find basic intersections (min of next group is smaller than max of current group & max of next group is larger than min of current group)
                    % In this instance, current group stops where next group stops
                    Tint1 = (Rmin(:,m+l) <= Rmax(:,m)).*(1-Tequ); Tint1 = (Tint1.*(Rmax(:,m+l) > Rmin(:,m))) > 0;
                    % Find advanced intersections (min/max of next group is smaller than min/max of current group) to deal with negative & 360 issues
                    % In this instance, current group starts where next group starts
                    Tint2 = Rmin(:,m+l) <= Rmin(:,m); Tint2 = (((Tint2 + (Rmax(:,m+l) <= Rmax(:,m))) > 0).*(1-Tequ)) > 0; 

                    % Check for & resolve intersects
                    if sum(Tint1)+sum(Tequ) == 3
                        % Update ranges
                        rx = (R(1,Bi_all(m))-ar(s)/2):a(s):(R(1,Bi_all(m+Tint1(1)*l))+ar(s)/2); % start at current group, stop at next group
                        ry = (R(2,Bi_all(m))-ar(s)/2):a(s):(R(2,Bi_all(m+Tint1(2)*l))+ar(s)/2);
                        rz = (R(3,Bi_all(m))-ar(s)/2):a(s):(R(3,Bi_all(m+Tint1(3)*l))+ar(s)/2);

                        % Add check
                        ccheck(m+l) = true;

                        % Update evaluators 
                        Rmax(:,m) = Rmax(:,m+l);
                        Rbase(:,m) = ((Rmax(:,m) - Rmin(:,m))./2) + Rmin(:,m);
                    elseif sum(Tint2)+sum(Tequ) == 3
                        % Correct for any wraps
                        corr = R(:,Bi_all(m+l)) - (Tint2.*360);

                        % Update ranges
                        rx = (corr(1)-ar(s)/2):a(s):(R(1,Bi_all(m))+ar(s)/2); % start at next group, start at current group
                        ry = (corr(2)-ar(s)/2):a(s):(R(2,Bi_all(m))+ar(s)/2);
                        rz = (corr(3)-ar(s)/2):a(s):(R(3,Bi_all(m))+ar(s)/2);

                        % Add check
                        ccheck(m+l) = true;

                        % Update evaluators
                        Rmin(:,m) = Rmin(:,m+l);
                        Rbase(:,m) = ((Rmax(:,m) - Rmin(:,m))./2) + Rmin(:,m);
                    end
                end
            end
        end

        % Obtain rotation angles
        Rx = repmat(rx, 1, length(ry)*length(rz)); Rx = sort(Rx);
        Ry = repmat(ry, 1, length(rz)); Ry = sort(Ry); Ry = repmat(Ry, 1, length(rx));
        Rz = repmat(rz, 1, length(rx)*length(ry));

        % Catalog new angles
        if m == 1
            Rnew = [Rx; Ry; Rz];
        else
            Rnew = horzcat(Rnew, [Rx; Ry; Rz]); %#ok<AGROW>
        end

        r(:,m) = [rx(1), ry(1), rz(1), rx(end), ry(end), rz(end)];
        r(r<0) = 360 + r(r<0);
        fprintf("      (%i, %i, %i) to (%i, %i, %i)\n", r(1,m), r(2,m), r(3,m), r(4,m), r(5,m), r(6,m))
    end

    % Correct for negative values & 360 values
    Rnew(Rnew<0) = 360 + Rnew(Rnew<0);
    Rnew(Rnew==360) = 0;

    % Add to old angles & remove duplicates
    if s == 1
        Rnew = DupAngleRmv(Rnew); R = Rnew;
    else
        Rold = R;
        R = DupAngleRmv(horzcat(R, Rnew));

        % Isolate new angles
        Rnew = R(:, (size(Rold,2)+1):end);
    end

    % Transform uCT data
    [m_test_uct] = Trans2Mask(i_test_uct, Rnew, Rin, ipmethod, thresh);

    % Update full-resolution check (log to see whether the angle has been evaluated at-full resolution)
    if size(m_test_uct,1) ~= size(m_fs_uct,1)
        fr = [fr; false(size(m_test_uct,3), 1)]; %#ok<AGROW>
    else
        fr = [fr; true(size(m_test_uct,3), 1)]; %#ok<AGROW>
    end

    % Evaluate rotatated FOV
    [dcetemp, ioutemp, RMSEtemp] = TransEval(m_test_uct, m_test_mri);
    a_iou = [a_iou; ioutemp]; a_dce = [a_dce; dcetemp]; a_RMSE = [a_RMSE; RMSEtemp]; %#ok<AGROW>

    % Find best FOV
    [~, Bi_dce] = maxk(a_dce, nume); % best dice values
    [~, Bi_iou] = maxk(a_iou, nume); % best iou values
    [~, Bi_RMSE] = mink(a_RMSE, nume); % best RMSE values

    % Get all "best" values
    Bi_all = union(Bi_dce, Bi_RMSE, 'stable');

    toc

    % Retransform & evaluate if any results listed as best were not found at full resolution
    if all(fr(Bi_all)) == false; fprintf('\nRe-testing the following angles at full resolution:\n'); end
    while all(fr(Bi_all)) == false
        % Find which angles were not tested at full resolution
        fri = intersect(find(fr == false), Bi_all);

        data = table(R(:,fri).', a_dce(fri), a_iou(fri), a_RMSE(fri), fr(fri), 'RowNames', string(fri));
        data.Properties.VariableNames = {'Rotation', 'dice', 'iou', 'RMSE', 'FReval'};
        data = sortrows(data, {'dice', 'RMSE'},{'descend', 'ascend'});
        disp(data)

        % Reset test mri data for re-evaluation
        m_test_mri = Norm01(sum(m_fs_mri,3));
        Rin = imref3d(size(i_fs_uct));

        % Re-transformation & re-evalution
        % Note: "R" is used instead of "Rnew" since indices listed in "B_all" are wrt "R"
        [m_test_uct] = Trans2Mask(i_fs_uct, R, Rin, 'cubic', thresh, fri);

        % Update full-resolution check
        fr(fri) = true(length(fri), 1);

        % Evaluate new data
        [dcetemp, ioutemp, RMSEtemp] = TransEval(m_test_uct, m_test_mri);
        a_iou(fri) = ioutemp;
        a_dce(fri) = dcetemp;
        a_RMSE(fri) = RMSEtemp;

        % Find best FOV
        [~, Bi_dce] = maxk(a_dce, nume); % best dice values
        [~, Bi_iou] = maxk(a_iou, nume); % best iou values
        [~, Bi_RMSE] = mink(a_RMSE, nume); % best RMSE values

        % Get all "best" values
        Bi_all = union(Bi_dce, Bi_RMSE, 'stable');
    end

    % Find overall best rotations
    Bi_all = intersect(Bi_dce, Bi_iou); Bi_all = intersect(Bi_all, Bi_RMSE);

    % Check for no intersections
    if isempty(Bi_all)
        % Create array of RMSE values for best dice/iou values
        RMSEtemp = ones(size(a_RMSE)); RMSEtemp(Bi_dce) = a_RMSE(Bi_dce);

        % Find max
        [~, Bi_all] = mink(RMSEtemp, numE(s));
    elseif length(Bi_all) > numE(s)
        Bi_all = Bi_all(1:numE(s));
    end

    data = table(R(:,Bi_all).', a_dce(Bi_all), a_iou(Bi_all), a_RMSE(Bi_all), 'RowNames', string(Bi_all));
    data.Properties.VariableNames = {'Rotation', 'dice', 'iou', 'RMSE'};
    fprintf('\nBest alignment occurred at rotation(s) showed below:\n')
    disp(data)

    toc
end

%% Show best FOV
% For display
B = Bi_all;
tform = affinetform3d(rigidtform3d(squeeze(R(:,B).'), [0 0 0]));

% Create MRI & uCT images for display
i_final_uct = imwarp(i_fs_uct, imref3d(size(i_fs_uct)), tform, 'cubic', "OutputView", affineOutputView(size(i_fs_uct), tform, "BoundsStyle","CenterOutput"));
if any(size(i_final_uct) ~= size(i_fs_uct)); i_temp = ROITrim(i_final_uct, size(i_fs_uct)); end
m_disp_uct = i_final_uct > thresh;
i_final_uct = i_final_uct.*m_disp_uct;
m_disp_uct = Norm01(sum(m_disp_uct, 3));
m_disp_uct = imfill(m_disp_uct);
i_disp_uct = Norm01(sum(i_final_uct, 3));
m_disp_mri = Norm01(sum(m_fs_mri, 3));
i_disp_mri = Norm01(sum(i_fs_mri, 3));

data = table(R.', a_dce, a_iou, a_RMSE, 'RowNames', string([1:length(a_iou)].')); %#ok<NBRAK1>
data.Properties.VariableNames = {'Rotation', 'dice', 'iou', 'RMSE'};

% take data off gpu for non-linear interpolation
if ~strcmp(ipmethod,'linear') && isgpuarray(i_fs_uct)
    i_fs_uct = gather(i_fs_uct);
end

% Show figures
figure('Name', 'mb', 'NumberTitle', 'off'); f1 = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
imagesc(m_disp_mri>0); axis square; axis off; colormap gray; title('MRI');
nexttile
imagesc(m_disp_uct>0); axis square; axis off; colormap gray; title('uCT');
nexttile
imagesc(imfuse(m_disp_mri>0, m_disp_uct>0)); axis square; axis off; title('Image Overlap'); hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'MRI Image','uCT Image'}); hold off;
title(f1, '2D Binary Masks')

figure('Name', 'm', 'NumberTitle', 'off'); f2 = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
imagesc(m_disp_mri); axis square; axis off; colormap gray; title('MRI');
nexttile
imagesc(m_disp_uct); axis square; axis off; colormap gray; title('uCT');
nexttile
imagesc(imfuse(m_disp_mri, m_disp_uct)); axis square; axis off; title('Image Overlap'); hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'MRI Image','uCT Image'}); hold off;
title(f2, '2D Masks')

figure('Name', 'i', 'NumberTitle', 'off'); f3 = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
imagesc(i_disp_mri); axis square; axis off; colormap gray; title('MRI');
nexttile
imagesc(i_disp_uct); axis square; axis off; colormap gray; title('uCT');
nexttile
imagesc(imfuse(i_disp_mri, i_disp_uct)); axis square; axis off; title('Image Overlap'); hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'MRI Image','uCT Image'}); hold off;
title(f3, '2D Images')

% dce_scale = a_dce.^4; dce_scale = (dce_scale./max(dce_scale)).*40;
% iou_scale = a_iou.^4; iou_scale = (iou_scale./max(iou_scale)).*40;
% RMSE_scale = a_RMSE.^4; RMSE_scale = (RMSE_scale./min(RMSE_scale)).*40;
% figure('Name', 'dice', 'NumberTitle', 'off'); scatter3(Rx,Ry,Rz,dce_scale, a_dce, 'filled'); xlim([0 360]); ylim([0 360]); zlim([0 360]); xlabel('X Rotation (Degrees)'); ylabel('Y Rotation (Degrees)'); zlabel('Z Rotation (Degrees)'); cb = colorbar; cb.Label.String = 'dice'; colormap(flipud(turbo));
% figure('Name', 'iou', 'NumberTitle', 'off'); scatter3(Rx,Ry,Rz,iou_scale, a_iou, 'filled'); xlim([0 360]); ylim([0 360]); zlim([0 360]); xlabel('X Rotation (Degrees)'); ylabel('Y Rotation (Degrees)'); zlabel('Z Rotation (Degrees)'); cb = colorbar; cb.Label.String = 'iou'; colormap(flipud(turbo));
% figure('Name', 'RMSE', 'NumberTitle', 'off'); scatter3(Rx,Ry,Rz,RMSE_scale, a_RMSE, 'filled'); xlim([0 360]); ylim([0 360]); zlim([0 360]); xlabel('X Rotation (Degrees)'); ylabel('Y Rotation (Degrees)'); zlabel('Z Rotation (Degrees)'); cb = colorbar; cb.Label.String = 'RMSE'; colormap(turbo);

% clear i j k i2 j2 m diff tissue CC idx numPixels i_temp m_temp tform im_name m_name m

%% Functions
% % Assuming both inputs are same FOV, decreases resolution of first image to match that of second image
% function img_newres = Downscale3D(img_oldres, imsize)
% dimlen = size(img_oldres); % find image size
% 
% % Check if input is square
% if min(dimlen) ~= max(dimlen)
% 
%     % Find size difference
%     diff = max(dimlen) - min(dimlen); diff1 = floor(diff/2);
% 
%     % Check if difference is odd (if so, remove extra from one side)
%     extra = rem(diff, 2);
% 
%     % Trim FOV
%     [~, Mdim] = max(dimlen);
%     if Mdim == 1
%         if dimlen(2) == dimlen(1)
%             img_oldres = img_oldres(1+diff1:dimlen-diff1-extra, 1+diff1:dimlen-diff1-extra, :);
%         elseif dimlen(3) == dimlen(1)
%             img_oldres = img_oldres(1+diff1:dimlen-diff1-extra, :, 1+diff1:dimlen-diff1-extra);
%         else
%             img_oldres = img_oldres(1+diff1:dimlen-diff1-extra, :, :);
%         end
%     elseif Mdim == 2
%         if dimlen(1) == dimlen(2)
%             img_oldres = img_oldres(1+diff1:dimlen-diff1-extra, 1+diff1:dimlen-diff1-extra, :);
%         elseif dimlen(3) == dimlen(2)
%             img_oldres = img_oldres(:, 1+diff1:dimlen-diff1-extra, 1+diff1:dimlen-diff1-extra);
%         else
%             img_oldres = img_oldres(:, 1+diff1:dimlen-diff1-extra, :);
%         end
%     else % Mdim ==
%         if dimlen(1) == dimlen(3)
%             img_oldres = img_oldres(1+diff1:dimlen-diff1-extra, :, 1+diff1:dimlen-diff1-extra);
%         elseif dimlen(2) == dimlen(3)
%             img_oldres = img_oldres(:, 1+diff1:dimlen-diff1-extra, 1+diff1:dimlen-diff1-extra);
%         else
%             img_oldres = img_oldres(:, :, 1+diff1:dimlen-diff1-extra);
%         end
%     end
% end
% 
% dimlen = size(img_oldres,1); % find image size
% 
% kspace = ifftshift(fftn(ifftshift(img_oldres))); % kspace of og data
% 
% % Find size difference
% diff = size(img_oldres, 1) - imsize(1); diff1 = floor(diff/2);
% 
% % Check if difference is odd (if so, remove extra from one side)
% extra = rem(diff, 2);
% 
% % Decrease resolution
% kspace = kspace(1+diff1:dimlen-diff1-extra, 1+diff1:dimlen-diff1-extra, 1+diff1:dimlen-diff1-extra);
% 
% % Convert from FD
% img_newres = abs(fftshift(ifftn(fftshift(kspace))));
% end
% 
% % Assuming both inputs are same FOV, decreases resolution of first image to match that of second image
% function img_newres = Downscale2D(img_oldres, imsize)
% dimlen = size(img_oldres,1); % find image size
% 
% kspace = ifftshift(fftn(ifftshift(img_oldres))); % kspace of og data
% 
% % Find size difference
% diff = size(img_oldres, 1) - imsize(1); diff1 = floor(diff/2);
% 
% % Check if difference is odd (if so, remove extra from one side)
% extra = rem(diff, 2);
% 
% % Decrease resolution
% kspace = kspace(1+diff1:dimlen-diff1-extra, 1+diff1:dimlen-diff1-extra, :);
% 
% % Convert from FD
% img_newres = abs(fftshift(ifftn(fftshift(kspace))));
% end

% Centers 2D & 3D data
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
X = (x - min(x,[],[1 2]))./max(x,[],[1 2]);
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

% q = 100;
% RMt = [cosd(R(2,q)).*cosd(R(3,q)), sind(R(1,q)).*sind(R(2,q)).*cosd(R(3,q))-cosd(R(1,q)).*sind(R(3,q)), cosd(R(1,q)).*sind(R(2,q)).*cosd(R(3,q))+sind(R(1,q)).*sind(R(3,q));
%        cosd(R(2,q)).*sind(R(3,q)), sind(R(1,q)).*sind(R(2,q)).*sind(R(3,q))+cosd(R(1,q)).*cosd(R(3,q)), cosd(R(1,q)).*sind(R(2,q)).*sind(R(3,q))-sind(R(1,q)).*cosd(R(3,q));
%        -sind(R(2,q)), sind(R(1,q)).*cosd(R(2,q)), cosd(R(1,q)).*cosd(R(2,q))];

% Correct for any rounding errors between GPU & CPU
RM = round(RM*10e9);

% Add first angle into eulcheck matrix
eulcheck(:,:,1) = RM(:,:,1);

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

% % Remove angles that would produce duplicate rotational matrices (4D data)
% function Rcorr = DupAngleRmv4D(R)
% % Check for unique rotational matrices
% Rcorr = zeros(size(R), 'like', R);
% eulcheck = zeros(size(R,1), size(R,1), size(R,2), size(R,3), size(R,4), 'like', R);
% RM = zeros(size(R,1), size(R,1), size(R,2), size(R,3), size(R,4), 'like', R);
% 
% RM(1,:,:,:,:) = [cosd(R(2,:,:,:)).*cosd(R(3,:,:,:)); sind(R(1,:,:,:)).*sind(R(2,:,:,:)).*cosd(R(3,:,:,:))-cosd(R(1,:,:,:)).*sind(R(3,:,:,:)); cosd(R(1,:,:,:)).*sind(R(2,:,:,:)).*cosd(R(3,:,:,:))+sind(R(1,:,:,:)).*sind(R(3,:,:,:))];
% RM(2,:,:,:,:) = [cosd(R(2,:,:,:)).*sind(R(3,:,:,:)); sind(R(1,:,:,:)).*sind(R(2,:,:,:)).*sind(R(3,:,:,:))+cosd(R(1,:,:,:)).*cosd(R(3,:,:,:)); cosd(R(1,:,:,:)).*sind(R(2,:,:,:)).*sind(R(3,:,:,:))-sind(R(1,:,:,:)).*cosd(R(3,:,:,:))];
% RM(3,:,:,:,:) = [-sind(R(2,:,:,:)); sind(R(1,:,:,:)).*cosd(R(2,:,:,:)); cosd(R(1,:,:,:)).*cosd(R(2,:,:,:))];
% 
% for i = 2:size(R,2)
%     for j = 1:size(R,3)
%         for k = 1:size(R,4)
%             % Check if rotational matrix is equal to any others already included
%             if all(sum(abs(eulcheck - RM(:,:,i,j,k)), [1 2]) ~= 0, 'all') && all(sum(abs(eulcheck + RM(:,:,i,j,k)), [1 2]) ~= 0, 'all')
%                 eulcheck(:,:,i,j,k) = RM(:,:,i,j,k);
%                 Rcorr(:,i,j,k) = R(:,i,j,k);
%             else
%                 Rcorr(:,i,j,k) = NaN(3,1);
%             end
%         end
%     end
% end
% end

% Transform image & convert to 2D mask
function [msk] = Trans2Mask(i_test, R, Rin, method, th, B)

% Clarify # iterations & iterator, based off input
switch nargin
    case 5
        inum = size(R,2); I = 1:inum;
    case 6
        inum = size(B,1); I = B;
    otherwise
        error('Incorrect number of input arguments.')
end

% Preallocate data
msk = zeros(size(i_test,1), size(i_test,2), inum, 'like', i_test);

for i = 1:inum
    % Create affine transform
    tform = affinetform3d(rigidtform3d(squeeze(R(:,I(i)).'), [0 0 0]));

    % Transform image
    % i_temp = imwarpcustom(i_test, Rin, tform, ipmethod, affineOutputView(size(i_test), tform, "BoundsStyle","CenterOutput"));
    i_temp = imwarp(i_test, Rin, tform, method, "OutputView", affineOutputView(size(i_test), tform, "BoundsStyle","CenterOutput"));

    % Trim to correct size (if applicable)
    if any(size(i_temp) ~= size(i_test)); i_temp = ROITrim(i_temp, size(i_test)); end

    % Create basic mask
    m_temp = i_temp > th;

    % Convert to 2D
    m_temp = sum(m_temp, 3);

    % Add into array
    msk(:,:,i) = m_temp;
end

% Normalize
msk = Norm01Slice(msk);
end

% Evaluate transformations
function [dce, iou, RMSE] = TransEval(m_tran, m_true)
m_true = repmat(m_true, [1 1 size(m_tran,3)]);

% RMSE
RMSE = sqrt(squeeze(sum((m_true - m_tran).^2,[1 2]))./(size(m_tran,1)*size(m_tran,2)));

% IOU & dice parameters
m_tran = m_tran > 0; m_true = m_true > 0;
area_i = squeeze(sum(m_tran.*m_true, [1 2]));
area_t = (squeeze(sum(m_tran, [1 2])) + squeeze(sum(m_true, [1 2])));

% IOU
iou = area_i./(area_t - area_i);
% dice
dce = 2.*area_i./area_t;
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
