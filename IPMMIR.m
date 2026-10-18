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
tissue = 'oregano';
res = 0.069447; % resolution, mm
thresh = 0.25; % uCT mask threshold
roipad = 20; % roi padding
ipmethod = 'linear';  % interpolation method
rsval = [  64,  96, 128     ]; % resize value
a =     [  30,  10,   2,   1,   1]; % rotation resolutions & ranges for transformation
ar =    [ 360,  60,  12,   6,   2];
numE =          [5,   4,   3,   1]; % number of local extrema from previous iteration to test
nume = 10; % number of local extrema before combination
F = length(a)-1; % number of iterations
W = 20; % max number of while loop iterations
dtype = 'image'; % type of data to use in test. Options are 'image' and 'mask'

%%% Load MRI data
[im_name, m_name, ~, se_offset] = grmdfs(tissue, 'ssems');
[~, mri_i_fs, P, ~, ~] = loadfidfancy(im_name);
% Import mask
load(m_name); mri_m_fs = m; 
% Circular shift (so it aligns w/ GRE)
mri_i_fs = circshift(mri_i_fs, se_offset, 1);
mri_m_fs = circshift(mri_m_fs, se_offset, 1);
% Apply mask
mri_i_fs = abs(mri_i_fs.*mri_m_fs);
% Translate to center
[diff, ~, ~] = CentShift(mri_i_fs, mri_m_fs);
mri_i_fs = circshift(mri_i_fs, [diff 0]); mri_m_fs = circshift(mri_m_fs, [diff 0]);
% Normalize data (done this way because 1,3,5,7 & 2,4,6,8,10 were acquired separately
mri_set1 = zeros(size(mri_i_fs,1), size(mri_i_fs,2), size(mri_i_fs,3)/2);
mri_set2 = zeros(size(mri_i_fs,1), size(mri_i_fs,2), size(mri_i_fs,3)/2);
I1 = 0; I2 = 0;
for i = 1:size(mri_i_fs,3) % split
    if rem(i,2) == 1
        I1 = I1 + 1;
        mri_set1(:,:,I1) = mri_i_fs(:,:,i);
    else
        I2 = I2 + 1;
        mri_set2(:,:,I2) = mri_i_fs(:,:,i);
    end
end
mri_set1 = Norm01(mri_set1); mri_set2 = Norm01(mri_set2); % normalize
I1 = 0; I2 = 0;
for i = 1:size(mri_i_fs,3) % recombine
    if rem(i,2) == 1
        I1 = I1 + 1;
        mri_i_fs(:,:,i) = mri_set1(:,:,I1);
    else
        I2 = I2 + 1;
        mri_i_fs(:,:,i) = mri_set2(:,:,I2);
    end
end

fov = P.FOV.fov.ro*10; % field of view, mm

%%% Import uCT data
uct_i_og = niftiread("C:\Users\jacob\OneDrive - Texas A&M University\TAMU\NMR RF Lab\_Project GRMD\uCT Data\SlicerSaves\Volumes\"+tissue+"_vol.nii");
% Convert data type to double
uct_i_og = double(uct_i_og);
% Reshape & rotate array to correct format (for some reason, it outputs flipped along both axes
% This moves the "coronal" slices to the sagittal position & rotates 180 along z-axis
uct_i_og = rot90(permute(uct_i_og, [3 1 2]),2); % 'up resolution' uCT image (in reality, just normal resolution)
% Ensure data is a cube of 23mm
for f1 = 1:3
    arrpad = [0 0 0];
    if size(uct_i_og, f1) < round(fov/res) % number of pixels for 23mm FOV of 69.447um resolution
        diff = round(fov/res) - size(uct_i_og, f1);
        if rem(diff,2) == 1
            diff = diff - 1; 
            arrpad(f1) = 1;
            uct_i_og = padarray(uct_i_og, arrpad, 0, 'pre');
        end
        arrpad(f1) = round(diff/2);
        uct_i_og = padarray(uct_i_og, arrpad, 0, 'both');
    end
end

resdiff = 1/(((res/size(mri_i_fs,1))*size(uct_i_og,1)));

clear diff f1 fov f2 m_name im_name arrpad P m se_offset

%% Create binary mask for uCT
% Downscale data
% uct_im_nt = downscale3D(uct_im_fr, size(mri_im_3D)); % 'not transformed' uCT image
uct_i_fs = Norm01(imresize3(uct_i_og, [size(mri_i_fs,1) size(mri_i_fs,1) size(mri_i_fs,1)]));

% Keep unmasked image for final analysis
uct_i_fsum = uct_i_fs;

% Create mask
uct_m_fs = GenMask(uct_i_fs, thresh);

% Apply to uCT image
uct_i_fs = uct_i_fsum.*uct_m_fs;

%% Trim ROI
mri_m_fs = ROITrim(mri_m_fs, roipad);
uct_m_fs = ROITrim(uct_m_fs, roipad);

% Check if same size (likely won't be)
if size(mri_m_fs,1) > size(uct_m_fs,1)
    roidiff = size(mri_m_fs,1) - size(uct_m_fs,1);

    % Pad array
    uct_m_fs = padarray(uct_m_fs, [round(roidiff/2) round(roidiff/2) round(roidiff/2)], 0, 'both');

    % Trim any excess
    if rem(roidiff, 2) > 0
        uct_m_fs = uct_m_fs(2:end,2:end,2:end);
    end
elseif size(mri_m_fs,1) < size(uct_m_fs,1)
    roidiff = size(uct_m_fs,1) - size(mri_m_fs,1);

    % Pad array
    mri_m_fs = padarray(mri_m_fs, [round(roidiff/2) round(roidiff/2) 0], 0, 'both');

    % Trim any excess
    if rem(roidiff, 2) > 0
        mri_m_fs = mri_m_fs(2:end,2:end,:);
    end
end

% Apply trim to images
mri_i_fs = ROITrim(mri_i_fs, size(mri_m_fs));
[uct_i_fs, roi] = ROITrim(uct_i_fs, size(uct_m_fs));

% Lastly, center data in 3D space
% Note: it is important to center the data after the ROI is trimmed. If
% done so beforehand, the ROI trimming will off-center the volume to ensure
% equal padding on all sides
[diff, ~, ~] = CentShift(squeeze(sum(mri_i_fs,3)), squeeze(sum(mri_m_fs,3))>0);
mri_i_fs = circshift(mri_i_fs, [diff 0]); mri_m_fs = circshift(mri_m_fs, [diff 0]);
[diff, uct_i_fs, uct_m_fs] = CentShift(uct_i_fs, uct_m_fs);

% Also apply shift & trim to unmasked data
uct_i_fsum = uct_i_fsum(roi(1,1):roi(1,2), roi(2,1):roi(2,2), roi(3,1):roi(3,2));
uct_i_fsum = circshift(uct_i_fsum, diff);

if strcmp(dtype, 'image')
    uct_test = uct_i_fs;
else
    uct_test = uct_m_fs;
end

fprintf('\nFinished loading & pre-processing data.\n')

toc

clear roidiff

%% Transform & Evaluate
% Initialize variables
MT_bi = 1; % best metric
RT = [180; 180; 180]; % rotations
MT = []; % evaluation metrics
FR = []; % full-resolution checks

for f1 = 1:F
    fprintf('\n--------------------------- Phase %i ---------------------------\n', f1)

    % Resize (if applicable)
    if f1 < F
        mri_m_test = Norm01(sum(imresize3(mri_m_fs, [rsval(f1) rsval(f1) size(mri_m_fs,3)]),3));
        uct_testin = Norm01(imresize3(uct_test, [rsval(f1) rsval(f1) rsval(f1)]));
    else
        mri_m_test = Norm01(sum(mri_m_fs,3));
        uct_testin = uct_test;
    end

    % Create spatial referencing information
    rin = imref3d(size(uct_testin));

    % Sort ranges according to index. Should make finding overlaps easier
    MT_bi = sort(MT_bi, 'ascend');

    fprintf('Testing the following degree range(s) with a %i degree step size:\n', a(f1))

    % Check for overlapping ranges
    ccheck = false(length(MT_bi), 1);

    % Find all rotation angles with no repeats
    for f2 = 1:length(MT_bi)
        % Skip overlapping groups already merged 
        if ccheck(f2) == true
            continue
        end

        % Obtain rotation angles
        RT_tmp1 = GenAngles(RT(:,MT_bi(f2)), a(f1), ar(f1));

        if f2 < length(MT_bi)
            for f3 = 1:(length(MT_bi)-f2)
                % Skip checks with already combined groups
                if ccheck(f2+f3) == true
                    continue
                else
                    % Obtain rotation angles for testing range
                    RT_tmp2 = GenAngles(RT(:,MT_bi(f2+f3)), a(f1), ar(f1));

                    % Check for intersection point
                    [i, ~, ~] = intersect(RT_tmp1.', RT_tmp2.', 'rows', 'stable');
                    if ~isempty(i)
                        % Add check
                        ccheck(f2+f3) = true;

                        % Update ranges (old technique adds significantly more angles to test)
                        RT_tmp1 = union(RT_tmp1.', RT_tmp2.', 'rows', 'stable');
                        RT_tmp1 = RT_tmp1.';
                    end
                end
            end
        end

        % Catalog new angles
        if f2 == 1
            RT_new = RT_tmp1;
        else
            RT_new = horzcat(RT_new, RT_tmp1); %#ok<AGROW>
        end

        RT_tmp1(RT_tmp1<1) = 360 + RT_tmp1(RT_tmp1<1);
        rt = [min(RT_tmp1(1,:)), min(RT_tmp1(2,:)), min(RT_tmp1(3,:)), max(RT_tmp1(1,:)), max(RT_tmp1(2,:)), max(RT_tmp1(3,:))];
        fprintf("      (%i, %i, %i) to (%i, %i, %i)\n", rt(1), rt(2), rt(3), rt(4), rt(5), rt(6))
    end

    % Add to old angles & remove duplicates
    if f1 == 1
        RT_new = DupAngleRmv(RT_new); RT = RT_new;
    else
        RT_old = RT;
        RT = DupAngleRmv(horzcat(RT, RT_new));

        % Isolate new angles
        RT_new = RT(:, (size(RT_old,2)+1):end);
    end

    fprintf("%i unique rotations will be performed.\n", size(RT_new,2))

    % Transform uCT data
    uct_m_testout = Trans2Mask(uct_testin, RT_new, rin, ipmethod, thresh);

    % Update full-resolution check (log to see whether the angle has been evaluated at-full resolution)
    FR = [FR; false(size(uct_m_testout,3), 1)]; %#ok<AGROW>

    % Evaluate rotatated FOV
    MT_tmp = custommetrics(uct_m_testout, uct_m_testout>0, mri_m_test, mri_m_test>0, 2);
    MT = [MT; MT_tmp]; %#ok<AGROW>

    % Find best FOV
    [~, MT_bi1] = maxk(MT(:,1), nume);
    [~, MT_bi2] = maxk(MT(:,2), nume);

    % Combine all best FOVs into one array
    MT_bi = union(MT_bi1, MT_bi2, 'stable'); 

    toc

    % Reset test mri data for re-evaluation
    mri_m_test = Norm01(sum(mri_m_fs,3));
    uct_testin = uct_test;
    rin = imref3d(size(uct_testin));

    % Retransform & evaluate if any results listed as best were not found at full resolution
    if all(FR(MT_bi)) == false; fprintf('\nRe-testing the following angles at full resolution:\n'); end
    while all(FR(MT_bi)) == false
        % Find which angles were not tested at full resolution
        nFRi = intersect(find(FR == false), MT_bi);

        data = table(RT(:,nFRi).', MT(nFRi,:), 'RowNames', string(nFRi));
        data.Properties.VariableNames = {'Rotation', 'Metrics'};
        data = sortrows(data, {'Metrics'},{'descend'});
        disp(data)

        % Re-transformation & re-evalution
        % Note: "RT" is used instead of "RT_new" since indices listed in "Bi_met" are wrt "R"
        uct_m_testout = Trans2Mask(uct_testin, RT, rin, 'cubic', thresh, nFRi);

        % Update full-resolution check
        FR(nFRi) = true(length(nFRi), 1);

        % Evaluate new data
        MT_tmp = custommetrics(uct_m_testout, uct_m_testout>0, mri_m_test, mri_m_test>0, 2);
        MT(nFRi,:) = MT_tmp;

        % Find best FOV
        [~, MT_bi1] = maxk(MT(:,1), nume);
        [~, MT_bi2] = maxk(MT(:,2), nume);

        % Combine all best FOVs into one array
        MT_bi = union(MT_bi1, MT_bi2, 'stable');
    end

    % Find best FOV
    MT_bi = IsolateBest(MT, RT, nume, numE(f1));

    data = table(RT(:,MT_bi).', MT(MT_bi,:), 'RowNames', string(MT_bi));
    data.Properties.VariableNames = {'Rotation', 'Metrics'};
    fprintf('\nMaxima occurred at rotation(s) showed below:\n')
    disp(data)
    
    toc
end

%% Final Phase

% Reset test mri data for final evaluation
mri_i_test = mri_i_fs;
mri_m_test = mri_m_fs;
uct_testin = uct_test;
rin = imref3d(size(uct_testin));
RT_rs = RT(:,MT_bi);
f2 = 1; % index
w = 1; % while loop index

% Calculate metric for resliced data of current rotation
[uct_i_testout, uct_m_testout] = RotateData(uct_testin, RT_rs, rin, thresh);
[uct_i_testout, uct_m_testout, ~, ~] = CT_reslicer_v2(uct_i_testout, uct_m_testout, mri_i_test, mri_m_test, resdiff);
MT_rs = custommetrics(uct_i_testout, uct_m_testout, mri_i_test, mri_m_test, 3);

% While loop that repeatedly checks the area around the "best" alignment until it stops moving
RT_b = [-5; -5; -5]; % set the maxima to something that would never occur
MT_bi_old = MT_bi; MT_bi = 1;
fprintf('\n--------------------------- Phase %i ---------------------------\n', f1+1)
fprintf('Testing the resliced volume around the maxima in a %i degree radius.\nThis will repeatedly occur until the maxima stops moving. Metrics @ current maxima:\n', ar(end)/2)
fprintf('%.4f @ (%i, %i, %i)\n', sum(MT_rs), RT_rs(1), RT_rs(2), RT_rs(3))
while any(RT_b ~= RT_rs(:,MT_bi)) && (w < W)
    RT_b = RT_rs(:,MT_bi);

    % Obtain rotation angles
    RT_new = GenAngles(RT_rs(:,MT_bi), a(end), ar(end));

    % Identify new angles
    if w == 1
        RT_old = RT_b;
    else
        RT_old = RT_rs;
    end
    RT_rs = horzcat(RT_old, RT_new);
    RT_rs = DupAngleRmv(RT_rs);

    % Isolate new angles
    RT_new = RT_rs(:, (size(RT_old,2)+1):end);

    % Preallocate data
    MT_rs = [MT_rs; zeros(size(RT_new,2),2)]; %#ok<AGROW>

    % Test new angles
    for f1 = 1:size(RT_new,2)
        f2 = f2 + 1;

        % Transform image & create mask
        [uct_i_testout, uct_m_testout] = RotateData(uct_testin, RT_new(:,f1), rin, thresh);

        % Reslice microCT data
        [uct_i_testout, uct_m_testout, ~, ~] = CT_reslicer_v2(uct_i_testout, uct_m_testout, mri_i_test, mri_m_test, resdiff);

        % Evaluate
        MT_rs(f2,:) = custommetrics(uct_i_testout, uct_m_testout, mri_i_test, mri_m_test, 3);
    end

    % Find best FOV
    MT_bi = IsolateBest(MT_rs, RT_rs, nume, numE(end));

    fprintf('%.4f @ (%i, %i, %i)\n', sum(MT_rs(MT_bi,:)), RT_rs(1,MT_bi), RT_rs(2,MT_bi), RT_rs(3,MT_bi))

    w = w + 1;
end

if w >= W
    fprintf('\nThe program could not converge within a reasonable time. Please make modifications and try again.\n')
end

%% Show best FOV
% For display
tform = rigidtform3d(squeeze(RT_rs(:,MT_bi).'), [0 0 0]);

% Create MRI & uCT images for display
[uct_i_FNum, uct_m_FN] = RotateData(uct_i_fsum, RT_rs(:,MT_bi), rin, thresh);
uct_i_FN = uct_i_FNum.*uct_m_FN;
[uct_i_FNrs, uct_m_FNrs, PIX_b, SLC_b] = CT_reslicer_v2(uct_i_FN, uct_m_FN, mri_i_fs, mri_m_fs, resdiff);
uct_m_tmp = Norm01(sum(uct_m_FN, 3));
uct_i_tmp = Norm01(sum(uct_i_FN, 3));
mri_m_tmp = Norm01(sum(mri_m_fs, 3));
mri_i_tmp = Norm01(sum(mri_i_fs, 3));

fprintf('\nFinal orientation: (%i, %i, %i)\nFinal slice offset: %i\nFinal pixel offset: %i\n', RT_rs(1,MT_bi), RT_rs(2,MT_bi), RT_rs(3,MT_bi), SLC_b, PIX_b);

toc

data = table(RT.', MT, FR, 'RowNames', string([1:length(MT)].')); %#ok<NBRAK1>
data.Properties.VariableNames = {'Rotation', 'Metrics', 'FReval'};

data_rs = table(RT_rs.', MT_rs, 'RowNames', string([1:length(MT_rs)].')); %#ok<NBRAK1>
data_rs.Properties.VariableNames = {'Rotation', 'Metrics'};

% take data off gpu for non-linear interpolation
if ~strcmp(ipmethod,'linear') && isgpuarray(uct_i_fs)
    uct_i_fs = gather(uct_i_fs);
end

% Show figures
figure('Name', 'mb', 'NumberTitle', 'off'); fig1 = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
imagesc(mri_m_tmp>0); axis square; axis off; colormap gray; title('MRI');
nexttile
imagesc(uct_m_tmp>0); axis square; axis off; colormap gray; title('uCT');
nexttile
imagesc(imfuse(mri_m_tmp>0, uct_m_tmp>0)); axis square; axis off; title('Image Overlap'); hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'MRI','uCT'}); hold off;
title(fig1, '2D Binary Masks')

figure('Name', 'm', 'NumberTitle', 'off'); fig2 = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
imagesc(mri_m_tmp); axis square; axis off; colormap gray; title('MRI');
nexttile
imagesc(uct_m_tmp); axis square; axis off; colormap gray; title('uCT');
nexttile
imagesc(imfuse(mri_m_tmp, uct_m_tmp)); axis square; axis off; title('Image Overlap'); hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'MRI','uCT'}); hold off;
title(fig2, '2D Masks')

figure('Name', 'i', 'NumberTitle', 'off'); fig3 = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
imagesc(mri_i_tmp); axis square; axis off; colormap gray; title('MRI');
nexttile
imagesc(uct_i_tmp); axis square; axis off; colormap gray; title('uCT');
nexttile
imagesc(imfuse(mri_i_tmp, uct_i_tmp)); axis square; axis off; title('Image Overlap'); hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'MRI','uCT'}); hold off;
title(fig3, '2D Images')

plotmygraph(mri_i_fs, 'MRI Images'); 

plotmygraph(uct_i_FNrs, 'microCT Images')

figure('Name', 'finali', 'NumberTitle', 'off'); fig4 = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
for f1 = 1:size(mri_i_fs,3)
    nexttile
    imagesc(imfuse(mri_i_fs(:,:,f1), uct_i_FNrs(:,:,f1))); axis square; axis off; title(num2str(f1));
end
hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'MRI','uCT'}); hold off;
title(fig4, "Image Slices, poff = "+num2str(PIX_b)+", soff = "+num2str(SLC_b))

figure('Name', 'finalm', 'NumberTitle', 'off'); fig5 = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
for f1 = 1:size(mri_m_fs,3)
    nexttile
    imagesc(imfuse(mri_m_fs(:,:,f1), uct_m_FNrs(:,:,f1))); axis square; axis off; title(num2str(f1));
end
hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'MRI','uCT'}); hold off;
title(fig5, "Mask Slices, poff = "+num2str(PIX_b)+", soff = "+num2str(SLC_b))

clear rt

%%
dce2D = dicecustom(uct_i_tmp, uct_m_tmp>0, mri_i_tmp, mri_m_tmp>0, 2);
dce3D = dicecustom(uct_i_FNrs, uct_m_FNrs, mri_i_fs, mri_m_fs, 3);
iou2D = ioucustom(uct_i_tmp, uct_m_tmp>0, mri_i_tmp, mri_m_tmp>0, 2);
iou3D = ioucustom(uct_i_FNrs, uct_m_FNrs, mri_i_fs, mri_m_fs, 3);
rmse2D = -1*irmsecustom(uct_i_tmp, uct_m_tmp>0, mri_i_tmp, mri_m_tmp>0, 2);
rmse3D = -1*irmsecustom(uct_i_FNrs, uct_m_FNrs, mri_i_fs, mri_m_fs, 3);

data = table(string(tissue), RT_rs(:,MT_bi).', SLC_b, PIX_b, dce2D, dce3D, iou2D, iou3D, rmse2D, rmse3D);
data.Properties.VariableNames = {'Tissue', 'Rotation', 'Slice Offset', 'Pixel Offset', 'DSE 2D', 'DSE 3D', 'IoU 2D', 'IoU 3D', 'RMSE 2D', 'RMSE 3D'};
disp(data)
%%
[optimizer, metric] = imregconfig("multimodal");
uct_i_FNrs_registered = imregister(uct_i_FNrs(:,:,4), mri_i_fs(:,:,4), 'rigid', optimizer, metric);
tform = imregtform(uct_i_FNrs(:,:,4), mri_i_fs(:,:,4), 'rigid', optimizer, metric);

figure; imagesc(imfuse(uct_i_FNrs(:,:,4), uct_i_FNrs_registered)); axis square; axis off; title('Registration'); hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'unreg','reg'}); hold off;
figure; imagesc(imfuse(mri_i_fs(:,:,4), uct_i_FNrs(:,:,4))); axis square; axis off; title('Registration'); hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'MRI','unreg'}); hold off;
figure; imagesc(imfuse(mri_i_fs(:,:,4), uct_i_FNrs_registered)); axis square; axis off; title('Registration'); hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'MRI','reg'}); hold off;

%% Functions

% Determine best result
function MT_bi = IsolateBest(met, RT, nume, numE)
% Find best
[~, MT_bi1] = maxk(met(:,1), nume);
[~, MT_bi2] = maxk(met(:,2), nume);

% Combine all best into one array
MT_bi = intersect(MT_bi1, MT_bi2, 'stable');

% Check if none meet requirements
if isempty(MT_bi)

    % Convert rotations to radians
    RT_b1 = deg2rad(RT(:,MT_bi1).');
    RT_b2 = deg2rad(RT(:,MT_bi2).');

    % Repeat comparison matrices to allow for every combination
    RT_b1rep = repmat(RT_b1, [size(RT_b1,1), 1]);
    RT_b2rep = repmat(RT_b2, [size(RT_b1,1), 1]);
    for k = 1:length(MT_bi1)-1
        RT_b1rep((k-1)*size(RT_b1,1)+1:k*size(RT_b1,1),:) = circshift(RT_b1, [k 0]);
    end

    % Compare rotations to identify which are closest
    RT_comparisons = sum(abs(angdiff(RT_b1rep, RT_b2rep)),2);

    % Isolate closest rotations
    [~, RT_closesti] = mink(RT_comparisons, numE);
    RT_closest = round(rad2deg(RT_b1rep(RT_closesti,:)));

    % Identify index of these rotations
    MT_bi = ismember(RT.', RT_closest, 'rows');
    MT_bi = find(MT_bi == true);
end

if length(MT_bi) > numE
    MT_bi = MT_bi(1:numE);
end
end

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
X = (x - min(x,[],'all'))./(max(x,[],'all') - min(x,[],'all'));
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
% First, correct for 0s, negatives & >360 values
R(R<1) = R(R<1) + 360;
R(R>360) = R(R>360) - 360;

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
Rcorr(:,1) = R(:,1);

% Check for duplicates
for i = 2:size(R,2)
    % Check if rotational matrix is equal to any others already included
    if all(sum(abs(eulcheck(:,:,1:i) - RM(:,:,i)), [1 2]) ~= 0, 'all') && all(sum(abs(eulcheck(:,:,1:i) + RM(:,:,i)), [1 2]) ~= 0, 'all') 
        I = I + 1;
        eulcheck(:,:,i) = RM(:,:,i);
        Rcorr(:,I) = R(:,i);
    end
end
Rcorr = Rcorr(:,1:I);
end

% Apply high resolution rotation to image data
function [IMG, MSK] = RotateData(img, R, Rin, th)

% Create affine transform
tform = rigidtform3d(squeeze(R.'), [0 0 0]);

% Transform image
% IMG = imwarpcustom(i_test, Rin, tform, ipmethod, affineOutputView(size(i_test), tform, "BoundsStyle","CenterOutput"));
IMG = imwarp(img, Rin, tform, 'cubic', "OutputView", affineOutputView(size(img), tform, "BoundsStyle", "CenterOutput"));

% Trim to correct size (if applicable)
if any(size(IMG) ~= size(img)); IMG = ROITrim(IMG, size(img)); end

% Normalize
IMG = Norm01(IMG);

% Create mask
MSK = GenMask(IMG, th);
end

% Transform image & convert to 2D mask
function dout = Trans2Mask(din, R, Rin, method, th, B)

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
dout = zeros(size(din,1), size(din,2), inum, 'like', din);

for i = 1:inum

    % Create affine transform
    tform = rigidtform3d(squeeze(R(:,I(i)).'), [0 0 0]);

    % Transform image
    % tmp = imwarpcustom(i_test, Rin, tform, ipmethod, affineOutputView(size(i_test), tform, "BoundsStyle","CenterOutput"));
    tmp = imwarp(din, Rin, tform, method, "OutputView", affineOutputView(size(din), tform, "BoundsStyle", "CenterOutput"));

    % Create mask
    if ~strcmp(method, 'cubic')
        tmp = tmp > th;
    else
        tmp = GenMask(tmp, th);
    end

    % Convert to 2D
    tmp = sum(tmp, 3);

    % Add into array
    dout(:,:,i) = tmp;
end

% Normalize
dout = (dout - min(dout,[],[1 2]))./(max(dout,[],[1 2]) - min(dout,[],[1 2]));
end

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
% tComp = tIntrinsictoWorldOutput/tform.A'*tWorldToIntrinsicInput;
% tComp(:,4)=[0 0 0 1]; % avoid round off issues due to inversion above
% 
% % Find the transform that takes from input intrinsic to output intrinsic
% tformComposite = affinetform3d(tComp');
% 
% % Form plaid grid of intrinsic points in output image.
% [dstXIntrinsic,dstYIntrinsic,dstZIntrinsic] = meshgrid(1:Rout.ImageSize(2), 1:Rout.ImageSize(1), 1:Rout.ImageSize(3));
% [srcXIntrinsic,srcYIntrinsic, srcZIntrinsic] = transformPointsForward(tformComposite, dstXIntrinsic, dstYIntrinsic, dstZIntrinsic);
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
