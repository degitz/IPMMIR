%% Import & Reassemble data
%%% Key
% rs: re-sized
% fs: full-sized
% og: original
%%%
%%%%%% Uses 2D binary & non-binary masks for translation analysis, 3D masks for slice analysis
%%%%%% Transforms image & finds mask, then converts mask to 2D. Outputs only mask
clear; clc; close all;

tic

% Parameters
tissue = 'lyle';
res = 0.069447; % resolution, mm
thresh = 0.25; % uCT mask threshold
roipad = 20; % roi padding
ipmethod = 'linear';  % interpolation method
rsval = [  64,  96, 128     ]; % resize value
a =     [  30,  10,   2,   1,   1]; % rotation resolutions & ranges for transformation
ar =    [ 360,  60,  12,   6,   2];
numE =  [        5,   4,   3]; % number of local extrema from previous iteration to test
S = length(a)-1; % number of iterations
etype = 'Dice';
etype_rs = etype;
W = 20; % max number of while loop iterations

% Add last index in numE. Just acts as placeholder as we only want 1 output angle
numE = [numE, 1];

% Select evaluation metric
switch etype
    case 'Dice'
        TransEval = @(v,w,x,y,z) dicecustom(v,w,x,y,z);
    case 'iRMSE'
        TransEval = @(v,w,x,y,z) irmsecustom(v,w,x,y,z);
    case 'IoU'
        TransEval = @(v,w,x,y,z) ioucustom(v,w,x,y,z);
    case'Dice_iRMSE'
        TransEval = @(v,w,x,y,z) dicecustom(v,w,x,y,z)+irmsecustom(v,w,x,y,z);
end
switch etype_rs
    case 'Dice'
        SliceEval = @(v,w,x,y,z) dicecustom(v,w,x,y,z);
    case 'iRMSE'
        SliceEval = @(v,w,x,y,z) irmsecustom(v,w,x,y,z);
    case 'IoU'
        SliceEval = @(v,w,x,y,z) ioucustom(v,w,x,y,z);
    case'Dice_iRMSE'
        SliceEval = @(v,w,x,y,z) dicecustom(v,w,x,y,z)+irmsecustom(v,w,x,y,z);
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
[diff, ~, ~] = CentShift(i_fs_mri, m_fs_mri);
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

clear diff i fov m m_name im_name arrpad P

%% Create binary mask for uCT
% Downscale data
% uct_im_nt = downscale3D(uct_im_fr, size(mri_im_3D)); % 'not transformed' uCT image
i_fs_uct = Norm01(imresize3(i_og_uct, [size(i_fs_mri,1) size(i_fs_mri,1) size(i_fs_mri,1)]));

% Keep unmasked image for final analysis
i_fsum_uct = i_fs_uct;

% Create mask
m_fs_uct = GenMask(i_fs_uct, thresh);

% Apply to uCT image
i_fs_uct = i_fsum_uct.*m_fs_uct;

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
[i_fs_uct, roi] = ROITrim(i_fs_uct, size(m_fs_uct));

% Lastly, center data in 3D space
% Note: it is important to center the data after the ROI is trimmed. If
% done so beforehand, the ROI trimming will off-center the volume to ensure
% equal padding on all sides
[diff, ~, ~] = CentShift(squeeze(sum(i_fs_mri,3)), squeeze(sum(m_fs_mri,3))>0);
i_fs_mri = circshift(i_fs_mri, [diff 0]); m_fs_mri = circshift(m_fs_mri, [diff 0]);
[diff, i_fs_uct, m_fs_uct] = CentShift(i_fs_uct, m_fs_uct);

% Also apply shift & trim to unmasked data
i_fsum_uct = i_fsum_uct(roi(1,1):roi(1,2), roi(2,1):roi(2,2), roi(3,1):roi(3,2));
i_fsum_uct = circshift(i_fsum_uct, diff);

fprintf('\nFinished loading & pre-processing data.\n')

toc

clear roidiff

%% Transform & Evaluate
Bi_met = 1; R = [180; 180; 180]; a_met = []; fr = []; % initial

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

    % Sort ranges according to index. Should make finding overlaps easier
    Bi_met = sort(Bi_met, 'ascend');

    fprintf('Testing the following degree range(s) with a %i degree step size:\n', a(s))

    % Check for overlapping ranges
    ccheck = false(length(Bi_met), 1);

    % Find all rotation angles with no repeats
    for m = 1:length(Bi_met)
        % Skip overlapping groups already merged 
        if ccheck(m) == true
            continue
        end

        % Obtain rotation angles
        Rtest1 = GenAngles(R(:,Bi_met(m)), a(s), ar(s));

        if m < length(Bi_met)
            for l = 1:(length(Bi_met)-m)
                % Skip checks with already combined groups
                if ccheck(m+l) == true
                    continue
                else
                    % Obtain rotation angles for testing range
                    Rtest2 = GenAngles(R(:,Bi_met(m+l)), a(s), ar(s));

                    % Check for intersection point
                    [Int, int1, int2] = intersect(Rtest1.', Rtest2.', 'rows', 'stable');
                    if ~isempty(Int)
                        % Add check
                        ccheck(m+l) = true;

                        % Update ranges (old technique adds significantly more angles to test)
                        Rtest1 = union(Rtest1.', Rtest2.', 'rows', 'stable');
                        Rtest1 = Rtest1.';
                    end
                end
            end
        end

        % Catalog new angles
        if m == 1
            Rnew = Rtest1;
        else
            Rnew = horzcat(Rnew, Rtest1); %#ok<AGROW>
        end

        r = [min(Rtest1(1,:)), min(Rtest1(2,:)), min(Rtest1(3,:)), max(Rtest1(1,:)), max(Rtest1(2,:)), max(Rtest1(3,:))];
        r(r<0) = 360 + r(r<0);
        fprintf("      (%i, %i, %i) to (%i, %i, %i)\n", r(1), r(2), r(3), r(4), r(5), r(6))
    end

    % Add to old angles & remove duplicates
    if s == 1
        Rnew = DupAngleRmv(Rnew); R = Rnew;
    else
        Rold = R;
        R = DupAngleRmv(horzcat(R, Rnew));

        % Isolate new angles
        Rnew = R(:, (size(Rold,2)+1):end);
    end

    fprintf("%i unique rotations will be performed.\n", size(Rnew,2))

    % Transform uCT data
    [m_test_uct] = Trans2Mask(i_test_uct, Rnew, Rin, ipmethod, thresh);

    % Update full-resolution check (log to see whether the angle has been evaluated at-full resolution)
    fr = [fr; false(size(m_test_uct,3), 1)]; %#ok<AGROW>

    % Evaluate rotatated FOV
    mettemp = TransEval(m_test_uct, m_test_uct>0, m_test_mri, m_test_mri>0, 2);
    a_met = [a_met; mettemp]; %#ok<AGROW>

    % Find best FOV
    [~, Bi_met] = maxk(a_met, numE(s));

    toc

    % Retransform & evaluate if any results listed as best were not found at full resolution
    if all(fr(Bi_met)) == false; fprintf('\nRe-testing the following angles at full resolution:\n'); end
    while all(fr(Bi_met)) == false
        % Find which angles were not tested at full resolution
        fri = intersect(find(fr == false), Bi_met);

        data = table(R(:,fri).', a_met(fri), 'RowNames', string(fri));
        data.Properties.VariableNames = {'Rotation', etype};
        data = sortrows(data, {etype},{'descend'});
        disp(data)

        % Reset test mri data for re-evaluation
        m_test_mri = Norm01(sum(m_fs_mri,3));
        Rin = imref3d(size(i_fs_uct));

        % Re-transformation & re-evalution
        % Note: "R" is used instead of "Rnew" since indices listed in "Bi_met" are wrt "R"
        [m_test_uct] = Trans2Mask(i_fs_uct, R, Rin, 'cubic', thresh, fri);

        % Update full-resolution check
        fr(fri) = true(length(fri), 1);

        % Evaluate new data
        mettemp = TransEval(m_test_uct, m_test_uct>0, m_test_mri, m_test_mri>0, 2);
        a_met(fri) = mettemp;

        % Find best FOV
        [~, Bi_met] = maxk(a_met, numE(s));
    end

    data = table(R(:,Bi_met).', a_met(Bi_met), 'RowNames', string(Bi_met));
    data.Properties.VariableNames = {'Rotation', etype};
    fprintf('\nMaxima occurred at rotation(s) showed below:\n')
    disp(data)

    toc
end

%% Final Phase

% Reset test mri data for final evaluation
B_R_old = R(:,Bi_met);
R2 = B_R_old;
Bi_met = 1;
m_test_mri = m_fs_mri; i_test_mri = i_fs_mri;
Rin = imref3d(size(i_fs_uct));
I = 1; % index
w = 0; % while loop index

% Calculate metric for resliced data of current rotation
tform = affinetform3d(rigidtform3d(squeeze(B_R_old.'), [0 0 0]));
% i_temp = imwarpcustom(i_fs_uct, Rin, tform, 'cubic', affineOutputView(size(i_fsum_uct), tform, "BoundsStyle","CenterOutput"));
i_temp = imwarp(i_fs_uct, Rin, tform, 'cubic', "OutputView", affineOutputView(size(i_fs_uct), tform, "BoundsStyle","CenterOutput"));
if any(size(i_temp) ~= size(i_fs_uct)); i_temp = ROITrim(i_temp, size(i_fs_uct)); end
i_temp = Norm01(i_temp);
m_temp = GenMask(i_temp, thresh);
[i_temp, m_temp, ~, ~] = CT_reslicer(i_temp, m_temp, i_test_mri, m_test_mri, 1/res, etype_rs);
mettemp = SliceEval(i_temp, m_temp, i_test_mri, m_test_mri, 3);
rs_met = mean(mettemp);

% While loop that repeatedly checks the area around the "best" alignment until it stops moving
B_R = [-5; -5; -5]; % set the maxima to something that would never occur
fprintf('\n--------------------------- Phase %i ---------------------------\n', s+1)
fprintf('Testing the resliced volume around the maxima in a %i degree radius.\nThis will repeatedly occur until the maxima stops moving. %s @ current maxima:\n', ar(end)/2, etype_rs)
fprintf('%.4f @ (%i, %i, %i)\n', rs_met, R2(1), R2(2), R2(3))
while any(B_R ~= R2(:,Bi_met)) && (w < W)
    B_R = R2(:,Bi_met);

    % Obtain rotation angles
    Rnew = GenAngles(R2(:,Bi_met), a(end), ar(end));

    % Identify new angles
    Rold = R2;
    R2 = DupAngleRmv(horzcat(R2, Rnew));
    Rnew = R2(:, (size(Rold,2)+1):end); % new angles

    % Preallocate data
    rs_met = [rs_met; zeros(size(Rnew,2),1)];  %#ok<AGROW>

    % Test new angles
    for i = 1:size(Rnew,2)
        I = I + 1;

        % Transform image & create mask
        tform = affinetform3d(rigidtform3d(squeeze(Rnew(:,i).'), [0 0 0]));
        % i_temp = imwarpcustom(i_fs_uct, Rin, tform, 'cubic', affineOutputView(size(i_test), tform, "BoundsStyle","CenterOutput"));
        i_temp = imwarp(i_fs_uct, Rin, tform, 'cubic', "OutputView", affineOutputView(size(i_fs_uct), tform, "BoundsStyle","CenterOutput"));
        if any(size(i_temp) ~= size(i_fs_uct)); i_temp = ROITrim(i_temp, size(i_fs_uct)); end
        i_temp = Norm01(i_temp);
        m_temp = GenMask(i_temp, thresh);

        % Reslice microCT data
        [i_temp, m_temp, ~, ~] = CT_reslicer(i_temp, m_temp, i_test_mri, m_test_mri, 1/res, etype_rs);

        % Evaluate
        rs_met(I) = SliceEval(i_temp, m_temp, i_test_mri, m_test_mri, 3);
    end

    % Find best FOV
    [~, Bi_met] = maxk(rs_met, numE(end));
    fprintf('%.4f @ (%i, %i, %i)\n', rs_met(Bi_met), R2(1,Bi_met), R2(2,Bi_met), R2(3,Bi_met))

    w = w + 1;
end

if w >= W
    fprintf('\nThe program could not converge within a reasonable time. Please make modifications and try again.\n')
end

fprintf('\nFinal orientation: (%i, %i, %i)\n', R2(1,Bi_met), R2(2,Bi_met), R2(3,Bi_met));

toc

%% Show best FOV
% For display
B = Bi_met;
tform = affinetform3d(rigidtform3d(squeeze(R2(:,B).'), [0 0 0]));

% Create MRI & uCT images for display
i_finalum_uct = imwarp(i_fsum_uct, imref3d(size(i_fsum_uct)), tform, 'cubic', "OutputView", affineOutputView(size(i_fsum_uct), tform, "BoundsStyle","CenterOutput"));
if any(size(i_finalum_uct) ~= size(i_fsum_uct)); i_temp = ROITrim(i_finalum_uct, size(i_fsum_uct)); end
i_finalum_uct = Norm01(i_finalum_uct);
m_final_uct = GenMask(i_finalum_uct, thresh);
i_final_uct = i_finalum_uct.*m_final_uct;
m_disp_uct = Norm01(sum(m_final_uct, 3));
i_disp_uct = Norm01(sum(i_final_uct, 3));
m_disp_mri = Norm01(sum(m_fs_mri, 3));
i_disp_mri = Norm01(sum(i_fs_mri, 3));

data = table(R.', a_met, fr, 'RowNames', string([1:length(a_met)].')); %#ok<NBRAK1>
data.Properties.VariableNames = {'Rotation', etype, 'FReval'};

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
imagesc(imfuse(m_disp_mri>0, m_disp_uct>0)); axis square; axis off; title('Image Overlap'); hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'MRI','uCT'}); hold off;
title(f1, '2D Binary Masks')

figure('Name', 'm', 'NumberTitle', 'off'); f2 = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
imagesc(m_disp_mri); axis square; axis off; colormap gray; title('MRI');
nexttile
imagesc(m_disp_uct); axis square; axis off; colormap gray; title('uCT');
nexttile
imagesc(imfuse(m_disp_mri, m_disp_uct)); axis square; axis off; title('Image Overlap'); hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'MRI','uCT'}); hold off;
title(f2, '2D Masks')

figure('Name', 'i', 'NumberTitle', 'off'); f3 = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
imagesc(i_disp_mri); axis square; axis off; colormap gray; title('MRI');
nexttile
imagesc(i_disp_uct); axis square; axis off; colormap gray; title('uCT');
nexttile
imagesc(imfuse(i_disp_mri, i_disp_uct)); axis square; axis off; title('Image Overlap'); hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'MRI','uCT'}); hold off;
title(f3, '2D Images')

%% Reslicing

% Reslice microCT data
[i_finaltrue_uct, m_finaltrue_uct, Bpix, Bofs] = CT_reslicer(i_final_uct, m_final_uct, i_fs_mri, m_fs_mri, 1/res, etype_rs);

plotmygraph(i_fs_mri, 'MRI Images'); plotmygraph(i_finaltrue_uct, 'microCT Images')
figure('Name', 'finali', 'NumberTitle', 'off'); f4 = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:size(i_fs_mri,3)
    nexttile
    imagesc(imfuse(i_fs_mri(:,:,i), i_finaltrue_uct(:,:,i))); axis square; axis off; title(num2str(i));
end
hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'MRI','uCT'}); hold off;
title(f4, "Image Slices, poff = "+num2str(Bpix)+", soff = "+num2str(Bofs))

plotmygraph(m_fs_mri, 'MRI Masks'); plotmygraph(m_finaltrue_uct, 'microCT Masks')
figure('Name', 'finalm', 'NumberTitle', 'off'); f5 = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:size(m_fs_mri,3)
    nexttile
    imagesc(imfuse(m_fs_mri(:,:,i), m_finaltrue_uct(:,:,i))); axis square; axis off; title(num2str(i));
end
hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'MRI','uCT'}); hold off;
title(f5, "Mask Slices, poff = "+num2str(Bpix)+", soff = "+num2str(Bofs))

%%
fprintf(['\n\ndice 2D: ', num2str(dicecustom(i_disp_uct, m_disp_uct>0, i_disp_mri, m_disp_mri>0, 2))])
fprintf(['\ndice 3D: ',num2str(dicecustom(i_finaltrue_uct, m_finaltrue_uct, i_fs_mri, m_fs_mri, 3))])
fprintf(['\nRMSE 2D: ',num2str(-1*irmsecustom(i_disp_uct, m_disp_uct>0, i_disp_mri, m_disp_mri>0, 2))])
fprintf(['\nRMSE 3D: ',num2str(-1*irmsecustom(i_finaltrue_uct, m_finaltrue_uct, i_fs_mri, m_fs_mri, 3)), '\n'])

%% Functions

% Centers 2D & 3D data
function [diff, imcs, varargout] = CentShift(im, m)
% Check for no mask input
if nargin == 1
   m = im;
   opcheck = false;
else
    opcheck = true;
end

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

% Check if we need to output mask as well
if opcheck
    varargout{1} = mcs;
end
end

% Normalize input data to range of 0-1
function X = Norm01(x)
X = (x - min(x,[],'all'))./max(x,[],'all');
end

% Generate binary mask
function m = GenMask(img, th)
% Create basic mask
m = img > th;

% Apply to uCT image
img = img.*m;

% Identify each connected region in image
CC = bwconncomp(m);

% Isolates largest region
[~,idx] = max(cellfun(@numel,CC.PixelIdxList)); 

% Create binary mask of only that region
m = false(size(img));
m(CC.PixelIdxList{idx}) = true;

% Fill any holes in mask
m = imfill(m, 6, "holes");
end

% Trim image down to specified size. Buff is either equal padding on all axis, or dimension size that roi should be padded to
function [T, roi] = ROITrim(ut, buff)
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

% Generate angle matrix from input ranges 
function R = GenAngles(r, a, ar)
rx = (r(1)-ar/2):a:(r(1)+ar/2);
ry = (r(2)-ar/2):a:(r(2)+ar/2);
rz = (r(3)-ar/2):a:(r(3)+ar/2);
Rx = repmat(rx, 1, length(ry)*length(rz)); Rx = sort(Rx);
Ry = repmat(ry, 1, length(rz)); Ry = sort(Ry); Ry = repmat(Ry, 1, length(rx));
Rz = repmat(rz, 1, length(rx)*length(ry));
R = [Rx; Ry; Rz];
end

% Remove angles that would produce duplicate rotational matrices
function Rcorr = DupAngleRmv(R)
% First, correct for negative values & +360 values
R(R<0) = R(R<0) + 360;
R(R>=360) = R(R>=360) - 360;

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
msk = (msk - min(msk,[],[1 2]))./max(msk,[],[1 2]);
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
% 
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
% 
% % Custom version of 'imwarp()'; more efficient & works better w/ GPU
% function outputImage = imwarpcustom(inputImage, Rin, tform, method, Rout)
% 
% % Define affine transformation that maps from intrinsic system of output image to world system of output image.
% tIntrinsictoWorldOutput = Rout.TransformIntrinsicToWorld;
% 
% % Define affine transformation that maps from world system of input image to intrinsic system of input image.
% tWorldToIntrinsicInput = Rin.TransformWorldToIntrinsic;
% 
% % Form transformation to go from output intrinsic to input intrinsic.
% % NOTE: tComp is in post-multiply form.
% tComp = tIntrinsictoWorldOutput/tform.A'*tWorldToIntrinsicInput; %tComp = Rout.TransformIntrinsicToWorld/tform.A'*Rin.TransformWorldToIntrinsic;
% tComp(:,4)=[0 0 0 1]; % avoid round off issues due to inversion above
% 
% % Find the transform that takes from input intrinsic to output intrinsic
% tformComposite = affinetform3d(tComp');
% 
% % Form plaid grid of intrinsic points in output image.
% [dstXIntrinsic,dstYIntrinsic,dstZIntrinsic] = meshgrid(1:Rout.ImageSize(2), 1:Rout.ImageSize(1), 1:Rout.ImageSize(3));
% [srcXIntrinsic,srcYIntrinsic, srcZIntrinsic] = transformPointsForward(tformComposite, dstXIntrinsic,dstYIntrinsic, dstZIntrinsic);
% 
% % Interpolate to produce transformed image
% outputImage = interp3(inputImage,srcXIntrinsic,srcYIntrinsic,srcZIntrinsic,method);
% 
% % Replace NaNs with zeros
% outputImage(isnan(outputImage)) = 0;
% 
% % outputImage = images.internal.interp3d(inputImage,srcXIntrinsic,srcYIntrinsic,srcZIntrinsic,method,fillValues, SmoothEdges);
% % works in a similar manner as interp3 when smoothing is turned off except it automatically replaces NaNs w/ zeros
% end
