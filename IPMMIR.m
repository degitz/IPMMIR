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
a      = [ 30,  10,   2,   1,   1]; % rotation resolutions for transformation
ar     = [360,  60,  12,   6,   2]; % rotation ranges for transformation
numE   =       [ 5,   4,   3,   1]; % number of local extrema from previous iteration to test
nume   = 10;                        % number of local extrema before combination
W      = 20;                        % max number of while loop iterations for final phase
roipad = 20;                        % roi padding, # pixels
th     =  0.25;                     % test mask threshold
dtype  = 'image';                   % type of data to use in test. Options are 'image' and 'mask'
tissue = 'dorothy';
F = 4;                    % number of iterations

%%% Load MRI data
[im_name, m_name, ~, se_offset] = grmdfs(tissue, 'ssems');
[~, mri_i_seum, P, ~, ~] = loadfidfancy('fid', im_name, 'zip', true);
load(m_name); mri_m_fs = m;
% Load GRE
[im_name, m_name, ~] = grmdfs(tissue, 'sgems');
[~, mri_i_geum, ~, ~, ~] = loadfidfancy('fid', im_name, 'zip', true); 
load(m_name); mri_m_ge = m; 
% Circular shift (so it aligns w/ GRE)
mri_i_seum = circshift(mri_i_seum, se_offset, 1);
mri_m_fs = circshift(mri_m_fs, se_offset, 1);
% Isolate magnitude
mri_i_seum = abs(mri_i_seum); mri_i_geum = abs(mri_i_geum); mri_i_geumog = mri_i_geum; mri_i_seumog = mri_i_seum;
% Translate to center
[diff_mri, ~, ~] = CentShift(mri_i_seum.*mri_m_fs, mri_m_fs);
mri_i_seum = circshift(mri_i_seum, [diff_mri 0]); mri_m_fs = circshift(mri_m_fs, [diff_mri 0]); mri_i_geum = circshift(mri_i_geum, [diff_mri 0]); mri_m_ge = circshift(mri_m_ge, [diff_mri 0]);
% Normalize data (done this way because 1,3,5,7 & 2,4,6,8,10 were acquired separately
mri_se_set1 = zeros(size(mri_i_seum,1), size(mri_i_seum,2), size(mri_i_seum,3)/2); mri_ge_set1 = mri_se_set1;
mri_se_set2 = zeros(size(mri_i_seum,1), size(mri_i_seum,2), size(mri_i_seum,3)/2); mri_ge_set2 = mri_se_set2;
I1 = 0; I2 = 0;
for i = 1:size(mri_i_seum,3) % split
    if rem(i,2) == 1
        I1 = I1 + 1;
        mri_se_set1(:,:,I1) = mri_i_seum(:,:,i); mri_ge_set1(:,:,I1) = mri_i_geum(:,:,i);
    else
        I2 = I2 + 1;
        mri_se_set2(:,:,I2) = mri_i_seum(:,:,i); mri_ge_set2(:,:,I2) = mri_i_geum(:,:,i);
    end
end
mri_se_set1 = Norm01(mri_se_set1); mri_se_set2 = Norm01(mri_se_set2); mri_ge_set1 = Norm01(mri_ge_set1); mri_ge_set2 = Norm01(mri_ge_set2); % normalize
I1 = 0; I2 = 0;
for i = 1:size(mri_i_seum,3) % recombine
    if rem(i,2) == 1
        I1 = I1 + 1;
        mri_i_seum(:,:,i) = mri_se_set1(:,:,I1); mri_i_geum(:,:,i) = mri_ge_set1(:,:,I1);
    else
        I2 = I2 + 1;
        mri_i_seum(:,:,i) = mri_se_set2(:,:,I2); mri_i_geum(:,:,i) = mri_ge_set2(:,:,I2);
    end
end
% Apply mask to final imported data
mri_i_fs = mri_i_seum.*mri_m_fs; mri_i_ge = mri_i_geum.*mri_m_ge;
% Isolate FOV size
fov = P.FOV.fov.ro*10; % field of view, mm

%%% Import uCT data
uct_i_og = niftiread("C:\Users\jacob\OneDrive - Texas A&M University\TAMU\NMR RF Lab\_Project GRMD\uCT Data\SlicerSaves\Volumes\"+tissue+"_vol.nii");
% Isolate resolution
res = niftiinfo("C:\Users\jacob\OneDrive - Texas A&M University\TAMU\NMR RF Lab\_Project GRMD\uCT Data\SlicerSaves\Volumes\"+tissue+"_vol.nii");
res = res.PixelDimensions(1); % resolution of test data, mm
% Convert data type to double
uct_i_og = double(uct_i_og);
% Reshape & rotate array to correct format (for some reason, it outputs flipped along both axes
% This moves the "coronal" slices to the sagittal position & rotates 180 along z-axis
uct_i_og = rot90(permute(uct_i_og, [3 1 2]),2); % 'up resolution' uCT image (in reality, just normal resolution)
% Ensure data is a cube of ref FOV
for f1 = 1:3
    arrpad = [0 0 0];
    if size(uct_i_og, f1) < round(fov/res) % number of pixels for ref FOV of test resolution
        diff_mri = round(fov/res) - size(uct_i_og, f1);
        if rem(diff_mri,2) == 1
            diff_mri = diff_mri - 1; 
            arrpad(f1) = 1;
            uct_i_og = padarray(uct_i_og, arrpad, 0, 'pre');
        end
        arrpad(f1) = round(diff_mri/2);
        uct_i_og = padarray(uct_i_og, arrpad, 0, 'both');
    end
end

clear diff_mri f1 fov f2 m_name im_name arrpad P m se_offset mri_se_set1 mri_se_set2 mri_ge_set1 mri_ge_set2

%% Determine test parameters
% Find closest power of 2
rsmax_uct = nextpow2(size(uct_i_og, 1)) - 1;
rsmax_mri = nextpow2(size(mri_i_fs, 1)) - 1;

% Keep whichever is smaller
if rsmax_mri < rsmax_uct
    rsmax = rsmax_mri;
else
    rsmax = rsmax_uct;
end

% Check if ref data is a power of 2
if (2^rsmax)*2 == size(mri_i_fs,1)
    rsmax = rsmax + 1;
end

% Calculate resize values
rsval = fliplr(round(2.^(rsmax-(1:F)./F.*2)));

clear rsmax_uct rsmax_mri rsmax

%% TEMPORARY FOR DATA CONSISTENCY
%%%%%%%%%%%%%%%%%%%%
rsval = [64, 96, 128, 256];
%%%%%%%%%%%%%%%%%%%%
%% Create binary mask for uCT
% Downscale data
% uct_im_nt = downscale3D(uct_im_fr, size(mri_im_3D)); % 'not transformed' uCT image
uct_i_fs = Norm01(imresize3(uct_i_og, [size(mri_i_fs,1) size(mri_i_fs,1) size(mri_i_fs,1)]));

% Calculate new resolution 
res_old = res;
res = res_old*(size(uct_i_og,1)/size(uct_i_fs,1));

% Difference in slice resolution (MRI slice thickness (mm) / uCT slice thickness (mm))
resdiff = 1/res;

% Create mask
uct_m_fs = GenMask(uct_i_fs, th);

% Keep unmasked image for final analysis
uct_i_um = uct_i_fs;

% Apply to uCT image
uct_i_fs = uct_i_um.*uct_m_fs;

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
[mri_i_fs, roi_mri] = ROITrim(mri_i_fs, size(mri_m_fs)); mri_i_ge = mri_i_ge(roi_mri(1,1):roi_mri(1,2), roi_mri(2,1):roi_mri(2,2), roi_mri(3,1):roi_mri(3,2)); mri_m_ge = mri_m_ge(roi_mri(1,1):roi_mri(1,2), roi_mri(2,1):roi_mri(2,2), roi_mri(3,1):roi_mri(3,2));
[uct_i_fs, roi_uct] = ROITrim(uct_i_fs, size(uct_m_fs));

% Lastly, center data in 3D space
% Note: it is important to center the data after the ROI is trimmed. If
% done so beforehand, the ROI trimming will off-center the volume to ensure
% equal padding on all sides
[diff_mri, ~, ~] = CentShift(squeeze(sum(mri_i_fs,3)), squeeze(sum(mri_m_fs,3))>0);
mri_i_fs = circshift(mri_i_fs, [diff_mri 0]); mri_m_fs = circshift(mri_m_fs, [diff_mri 0]); mri_i_ge = circshift(mri_i_ge, [diff_mri 0]); mri_m_ge = circshift(mri_m_ge, [diff_mri 0]);
[diff_uct, uct_i_fs, uct_m_fs] = CentShift(uct_i_fs, uct_m_fs);

% Also apply shift & trim to unmasked data
mri_i_seum = mri_i_seum(roi_mri(1,1):roi_mri(1,2), roi_mri(2,1):roi_mri(2,2), roi_mri(3,1):roi_mri(3,2)); mri_i_seum = circshift(mri_i_seum, [diff_mri 0]);
mri_i_geum = mri_i_geum(roi_mri(1,1):roi_mri(1,2), roi_mri(2,1):roi_mri(2,2), roi_mri(3,1):roi_mri(3,2)); mri_i_geum = circshift(mri_i_geum, [diff_mri 0]);
uct_i_um = uct_i_um(roi_uct(1,1):roi_uct(1,2), roi_uct(2,1):roi_uct(2,2), roi_uct(3,1):roi_uct(3,2)); uct_i_um = circshift(uct_i_um, diff_uct);

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
MT_bi = 1;            % best metric
RT = [180; 180; 180]; % rotations
MT = [];              % evaluation metrics
FR = [];              % full-resolution checks

for f1 = 1:F
    fprintf('\n--------------------------- Phase %i ---------------------------\n', f1)

    % Resize (if applicable)
    if size(mri_m_fs,1) ~= rsval(f1)
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
    uct_m_testout = Trans2Mask(uct_testin, RT_new, rin, 'linear', th);

    % Update full-resolution check (log to see whether the angle has been evaluated at-full resolution)
    FR = [FR; false(size(uct_m_testout,3), 1)]; %#ok<AGROW>

    % Evaluate rotatated FOV
    MT_tmp = custommetrics(uct_m_testout, uct_m_testout>0, mri_m_test, mri_m_test>0, 2);
    MT = [MT; MT_tmp]; %#ok<AGROW>

    % Find best FOV
    [~, MT_bi1] = maxk(MT(:,1), nume);
    [~, MT_bi2] = maxk(MT(:,2), nume);

    % Remove close extrema
    MT_bi1 = CloseExtrRmv(MT_bi1, RT, a(f1)*2);
    MT_bi2 = CloseExtrRmv(MT_bi2, RT, a(f1)*2);

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
        uct_m_testout = Trans2Mask(uct_testin, RT, rin, 'cubic', th, nFRi);

        % Update full-resolution check
        FR(nFRi) = true(length(nFRi), 1);

        % Evaluate new data
        MT_tmp = custommetrics(uct_m_testout, uct_m_testout>0, mri_m_test, mri_m_test>0, 2);
        MT(nFRi,:) = MT_tmp;

        % Reset best metrics indices
        MT_bi = [];

        % Create temporary extrema number
        nume_tmp = nume/1.5;

        % Ensure enough extrema are isolated
        % Note: DO NOT REMOVE OR MODIFY THE SECOND CONDITION. It is in 
        % place to ensure this doesn't become an infinite loop.
        while (length(MT_bi) < nume) && (nume_tmp < size(MT,1))
            % Increase number of extrema isolated
            nume_tmp = round(nume_tmp*1.5);

            % Find best FOV
            [~, MT_bi1] = maxk(MT(:,1), nume_tmp);
            [~, MT_bi2] = maxk(MT(:,2), nume_tmp);

            % Remove close extrema
            MT_bi1 = CloseExtrRmv(MT_bi1, RT, a(f1)*2);
            MT_bi2 = CloseExtrRmv(MT_bi2, RT, a(f1)*2);

            % Combine all best FOVs into one array
            MT_bi = union(MT_bi1, MT_bi2, 'stable');
        end
    end

    % Find best FOV
    MT_bi = IsolateBest(MT, RT, nume, numE(f1));

    data = table(RT(:,MT_bi).', MT(MT_bi,:), 'RowNames', string(MT_bi));
    data.Properties.VariableNames = {'Rotation', 'Metrics'};
    fprintf('\nMaxima occurred at rotation(s) showed below:\n')
    disp(data)
    
    toc
end

MT_bi_old = MT_bi;

%% Final Phase

% Reset test mri data for final evaluation
mri_i_test = mri_i_fs;
mri_m_test = mri_m_fs;
uct_testin = uct_test;
rin = imref3d(size(uct_testin));
RT_rs = RT(:,MT_bi_old);
f2 = 1; % index
w = 1; % while loop index

% Calculate metric for resliced data of current rotation
[uct_i_testout, uct_m_testout] = RotateData(uct_testin, RT_rs, rin, th);
[uct_i_testout, uct_m_testout, ~, ~] = Reslice(uct_i_testout, uct_m_testout, resdiff, mri_i_test, mri_m_test);
MT_rs = MII(uct_i_testout, mri_i_test);

% While loop that repeatedly checks the area around the "best" alignment until it stops moving
RT_b = [-5; -5; -5]; % set the maxima to something that would never occur
MT_bi_old = MT_bi; MT_bi = 1;
fprintf('\n--------------------------- Phase %i ---------------------------\n', f1+1)
fprintf('Testing the resliced volume around the maxima in a %i degree radius.\nThis will repeatedly occur until the maxima stops moving. Metrics @ current maxima:\n', ar(end)/2)
fprintf('%.4f @ (%i, %i, %i)\n', MT_rs, RT_rs(1), RT_rs(2), RT_rs(3))
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
    MT_rs = [MT_rs; zeros(size(RT_new,2),1)]; %#ok<AGROW>

    % Test new angles
    for f1 = 1:size(RT_new,2)
        f2 = f2 + 1;

        % Transform image & create mask
        [uct_i_testout, uct_m_testout] = RotateData(uct_testin, RT_new(:,f1), rin, th);

        % Reslice microCT data
        [uct_i_testout, uct_m_testout, ~, ~] = Reslice(uct_i_testout, uct_m_testout, resdiff, mri_i_test, mri_m_test);

        % Evaluate
        MT_rs(f2) = MII(uct_i_testout, mri_i_test);
    end

    % Find best FOV
    [~, MT_bi] = maxk(MT_rs, numE(end));

    fprintf('%.4f @ (%i, %i, %i)\n', sum(MT_rs(MT_bi)), RT_rs(1,MT_bi), RT_rs(2,MT_bi), RT_rs(3,MT_bi))

    w = w + 1;
end

if w >= W
    fprintf('\nThe program could not converge within a reasonable time. Please make modifications and try again.\n')
end

%% Show best FOV
% Create MRI & uCT images for display
[uct_i_FNum, uct_m_FN] = RotateData(uct_i_um, RT_rs(:,MT_bi), rin, th);
uct_i_FN = uct_i_FNum.*uct_m_FN;
[uct_i_FNrs, uct_m_FNrs, PIX_b, SLC_b] = Reslice(uct_i_FN, uct_m_FN, resdiff, mri_i_fs, mri_m_fs);
uct_m_tmp = Norm01(sum(uct_m_FN, 3));
uct_i_tmp = Norm01(sum(uct_i_FN, 3));
mri_m_tmp = Norm01(sum(mri_m_fs, 3));
mri_i_tmp = Norm01(sum(mri_i_fs, 3));

data = table(RT.', MT, FR, 'RowNames', string([1:length(MT)].')); %#ok<NBRAK1>
data.Properties.VariableNames = {'Rotation', 'Metrics', 'FReval'};

data_rs = table(RT_rs.', MT_rs, 'RowNames', string([1:length(MT_rs)].')); %#ok<NBRAK1>
data_rs.Properties.VariableNames = {'Rotation', 'Metrics'};

% take data off gpu for non-linear interpolation
if isgpuarray(uct_i_fs)
    uct_i_fs = gather(uct_i_fs);
end

plotmygraph(mri_i_seum, 'MRI SE UM Images'); 

plotmygraph(mri_i_geum, 'MRI GRE UM Images'); 

plotmygraph(mri_i_fs, 'MRI SE Images'); 

plotmygraph(mri_i_ge, 'MRI GRE Images'); 

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

% Output final metrics
dce2D = dicecustom(uct_i_tmp, uct_m_tmp>0, mri_i_tmp, mri_m_tmp>0, 2);
dce3D = dicecustom(uct_i_FNrs, uct_m_FNrs, mri_i_fs, mri_m_fs, 3);
iou2D = ioucustom(uct_i_tmp, uct_m_tmp>0, mri_i_tmp, mri_m_tmp>0, 2);
iou3D = ioucustom(uct_i_FNrs, uct_m_FNrs, mri_i_fs, mri_m_fs, 3);
rmse2D = -1*irmsecustom(uct_i_tmp, uct_m_tmp>0, mri_i_tmp, mri_m_tmp>0, 2);
rmse3D = -1*irmsecustom(uct_i_FNrs, uct_m_FNrs, mri_i_fs, mri_m_fs, 3);

run_metrics = table(sum(FR), size(RT,2), size(RT_rs,2));
run_metrics.Properties.VariableNames = {'Full Resolution Evaluations', 'Base Rotations', 'Resliced Rotations'};
disp(run_metrics)

final_metrics = table(string(tissue), RT_rs(:,MT_bi).', SLC_b, PIX_b, dce2D, dce3D, iou2D, iou3D, rmse2D, rmse3D);
final_metrics.Properties.VariableNames = {'Tissue', 'Rotation', 'Slice Offset', 'Pixel Offset', 'DSE 2D', 'DSE 3D', 'IoU 2D', 'IoU 3D', 'RMSE 2D', 'RMSE 3D'};
disp(final_metrics)

toc

% %%
% [optimizer, metric] = imregconfig("multimodal");
% uct_i_FNrs_registered = imregister(uct_i_FNrs(:,:,4), mri_i_fs(:,:,4), 'rigid', optimizer, metric);
% tform = imregtform(uct_i_FNrs(:,:,4), mri_i_fs(:,:,4), 'rigid', optimizer, metric);
% 
% figure; imagesc(imfuse(uct_i_FNrs(:,:,4), uct_i_FNrs_registered)); axis square; axis off; title('Registration'); hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'unreg','reg'}); hold off;
% figure; imagesc(imfuse(mri_i_fs(:,:,4), uct_i_FNrs(:,:,4))); axis square; axis off; title('Registration'); hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'MRI','unreg'}); hold off;
% figure; imagesc(imfuse(mri_i_fs(:,:,4), uct_i_FNrs_registered)); axis square; axis off; title('Registration'); hold on; legend([scatter([],[],10,[0 1 0],'filled') scatter([],[],10,[1 0 1],'filled')],{'MRI','reg'}); hold off;
% 
% I = mutualInfo(uct_i_FNrs(:,:,4), mri_i_fs(:,:,4));

%% Functions

% Centers 2D & 3D data
function [diff, imcs, varargout] = CentShift(im, m)
% Inputs
%       im: image
%        m: mask
% Outputs
%     diff: shift amount
%     imcs: centered image
%      mcs: shifted mask

% Check for no mask input
if nargin == 1
   m = im;
   flagM = false;
else
    flagM = true;
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
if flagM
    varargout{1} = mcs;
end
end

% Normalize input data to range of 0-1
function X = Norm01(x)
% Inputs
%        x: un-normalized data
% Outputs
%        X: normalized data

X = (x - min(x,[],'all'))./(max(x,[],'all') - min(x,[],'all'));
end

% Generate binary mask
function m = GenMask(img, th)
% Inputs
%      img: image
%       th: binary threshold
% Outputs
%        m: mask

% Create basic mask
m = img > th;

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

% Trim image down to specified size
function [T, roi] = ROITrim(ut, buff)
% Inputs
%       ut: un-trimmed image data
%     buff: either equal padding on all axis, or dimension size that roi should be padded to
% Outputs
%        T: shift amount
%      roi: trimmed image data

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
% Inputs
%       r: initial rotation
%       a: rotation range resolution
%      ar: rotation range
% Outputs
%       R: rotations

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
% Inputs
%        R: rotation matrix
% Outputs
%    Rcorr: corrected rotation matrix

% First, correct for 0s, negatives & >360 values
R(R<1) = R(R<1) + 360;
R(R>360) = R(R>360) - 360;

Rcorr = zeros(size(R), 'like', R); I = 1;
eulcheck = zeros(size(R,1), size(R,1), size(R,2), 'like', R);
RM = zeros(size(R,1), size(R,1), size(R,2), 'like', R);

% Calculate rotational matrices
RM(1,:,:) = [cosd(R(2,:)).*cosd(R(3,:)); sind(R(1,:)).*sind(R(2,:)).*cosd(R(3,:))-cosd(R(1,:)).*sind(R(3,:)); cosd(R(1,:)).*sind(R(2,:)).*cosd(R(3,:))+sind(R(1,:)).*sind(R(3,:))];
RM(2,:,:) = [cosd(R(2,:)).*sind(R(3,:)); sind(R(1,:)).*sind(R(2,:)).*sind(R(3,:))+cosd(R(1,:)).*cosd(R(3,:)); cosd(R(1,:)).*sind(R(2,:)).*sind(R(3,:))-sind(R(1,:)).*cosd(R(3,:))];
RM(3,:,:) = [-sind(R(2,:));              sind(R(1,:)).*cosd(R(2,:));                                          cosd(R(1,:)).*cosd(R(2,:))];

% q = 100;
% RMt = [cosd(R(2,q)).*cosd(R(3,q)), sind(R(1,q)).*sind(R(2,q)).*cosd(R(3,q))-cosd(R(1,q)).*sind(R(3,q)), cosd(R(1,q)).*sind(R(2,q)).*cosd(R(3,q))+sind(R(1,q)).*sind(R(3,q));
%        cosd(R(2,q)).*sind(R(3,q)), sind(R(1,q)).*sind(R(2,q)).*sind(R(3,q))+cosd(R(1,q)).*cosd(R(3,q)), cosd(R(1,q)).*sind(R(2,q)).*sind(R(3,q))-sind(R(1,q)).*cosd(R(3,q));
%        -sind(R(2,q)), sind(R(1,q)).*cosd(R(2,q)), cosd(R(1,q)).*cosd(R(2,q))];
% 
% R = RT_rs(:,q).';
% Rx = [1 0 0; 0 cosd(R(1)) -sind(R(1)); 0 sind(R(1)) cosd(R(1))];
% Ry = [cosd(R(2)) 0 sind(R(2)); 0 1 0; -sind(R(2)) 0 cosd(R(2))];
% Rz = [cosd(R(3)) -sind(R(3)) 0; sind(R(3)) cosd(R(3)) 0; 0 0 1];
% RMt2 = Rz*Ry*Rx;

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

% Remove extrema that are too close
function bi = CloseExtrRmv(bi, RT, th)
% Inputs
%       bi: best indices
%       RT: entire rotation matrix
%       th: threshold
% Outputs
%        bi: best indices, with close rotations removed

% Initialize increment
i1 = 1;

% Begin evaluation
while i1 < length(bi)
    % Initialize first quaternion
    q1 = quaternion(RT(:,bi(i1)).', 'eulerd', 'XYZ', 'point');

    % Compare distance from other extrema
    for i2 = (i1+1):length(bi)
        % Initialize second quaternion
        q2 = quaternion(RT(:,bi(i2)).', 'eulerd', 'XYZ', 'point');

        % Calculate distance
        qdist = dist(q1, q2)*180/pi;

        % Replace with NaN if distance is less than threshold
        if qdist < th
            bi(i2) = NaN;
        end
    end

    % Remove close extrema
    bi = rmmissing(bi);

    % Shift increment
    i1 = i1 + 1;
end
end

% Determine best result
function MT_bi = IsolateBest(MT, RT, nume, numE)
% Inputs
%       MT: metrixs
%       RT: rotations
%     nume: original number of local extrema
%     numE: final number of local extrema
% Outputs
%       bi: indices corresponding to the best metrics

% Find best (rotations which have highest DICE & RMSE)
[~, MT_bi1] = maxk(MT(:,1), nume);
[~, MT_bi2] = maxk(MT(:,2), nume);

% Isolate common rotations between the two arrays
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

% Transform image & convert to 2D mask
function MSK = Trans2Mask(img, R, Rin, method, th, idx)
% Inputs
%      img: image data
%        R: rotation(s)
%      Rin: spatial referencing information
%   method: interpolation method
%       th: mask threshold
%      idx: *optional* exact indices from rotation matrix to use
% Outputs
%      MSK: warped mask

% Clarify # iterations & iterator, based off input
switch nargin
    case 5
        inum = size(R,2); I = 1:inum;
    case 6
        inum = size(idx,1); I = idx;
    otherwise
        error('Incorrect number of input arguments.')
end

% Preallocate data
MSK = zeros(size(img,1), size(img,2), inum, 'like', img);

for i = 1:inum

    % Create affine transform
    tform = rigidtform3d(squeeze(R(:,I(i)).'), [0 0 0]);

    % Transform image
    % tmp = imwarpcustom(i_test, Rin, tform, ipmethod, affineOutputView(size(i_test), tform, "BoundsStyle","CenterOutput"));
    tmp = imwarp(img, Rin, tform, method, "OutputView", affineOutputView(size(img), tform, "BoundsStyle", "CenterOutput"));

    % Create mask
    if ~strcmp(method, 'cubic')
        tmp = tmp > th;
    else
        tmp = GenMask(tmp, th);
    end

    % Convert to 2D
    tmp = sum(tmp, 3);

    % Add into array
    MSK(:,:,i) = tmp;
end

% Normalize
MSK = (MSK - min(MSK,[],[1 2]))./(max(MSK,[],[1 2]) - min(MSK,[],[1 2]));
end

% Apply high resolution rotation to image data
function [IMG, MSK] = RotateData(img, R, Rin, th)
% Inputs
%      img: image data
%        R: rotation(s)
%      Rin: spatial referencing information
%       th: mask threshold
% Outputs
%      IMG: warped image
%      MSK: warped mask

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

% Reslice uCT data to match MRI slices
function [IMG, MSK, varargout] = Reslice(img, msk, r, inp4, inp5)
% Inputs
%      img: warped image data
%      msk: warped mask data
%        r: ratio of slice resolution between reference & test data
%     inp4: either pixel offset (prior data), or reference image data
%     inp5: either slice offset (prior data), or reference mask data
% Outputs
%      IMG: warped, resliced image
%      MSK: warped, resliced mask
%     Bpix: *optional* best pixel offset
%     Bslc: *optional* best slice offset

% Determine indices of start/stop for mri slices in ct
ibase = 0:r:size(img,3); % base starts & stops
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
    IMG = SliceCombine(img, ifull+Bpix, ipart);
    MSK = SliceCombine(msk, ifull+Bpix, ipart);
    IMG = IMG(:,:,firstslice:lastslice);
    MSK = MSK(:,:,firstslice:lastslice);
else
    img_ref = inp4; msk_ref = inp5;

    % Determine size difference of data
    slcdiff = size(ifull,2)-size(img_ref,3);
    pixdiff = size(img,3)-ifull(2,end);

    % Preallocate data
    MTpix = zeros(pixdiff+1, 2);
    Bslc = zeros(pixdiff+1, 1);

    % Compare CT slice assignments with available subslice offsets. Last offset
    % is excluded, since partial pixel addition would extend past matrix bounds.
    for pixofs = 0:pixdiff
        % fprintf('\nTesting a subslice offset of %i pixel(s).\n', pix)
        IMG = SliceCombine(img, ifull+pixofs, ipart);
        MSK = SliceCombine(msk, ifull+pixofs, ipart);

        % For extra slices in CT data, determine offset with MRI
        % Preallocate data
        MTslc = zeros(slcdiff+1, 2);
        for slcofs = 0:slcdiff

            % Isolate current number of slices
            i_CT_ofs = IMG(:,:,1+slcofs:end-slcdiff+slcofs);
            m_CT_ofs = MSK(:,:,1+slcofs:end-slcdiff+slcofs);

            MTslc(1+slcofs,:) = custommetrics(i_CT_ofs, m_CT_ofs, img_ref, msk_ref, 3);
        end

        % Find best
        Bslc(1+pixofs) = IsolateBestRS(MTslc)-1;

        % Isolate best from the rest
        IMG = IMG(:,:,1+Bslc(1+pixofs):end-slcdiff+Bslc(1+pixofs));
        MSK = MSK(:,:,1+Bslc(1+pixofs):end-slcdiff+Bslc(1+pixofs));

        % Evaluate
        MTpix(1+pixofs,:) = custommetrics(IMG, MSK, img_ref, msk_ref, 3);
    end

    % Find best
    Bpix = IsolateBestRS(MTpix);

    % Isolate best from the rest
    Bslc = Bslc(Bpix); Bpix = Bpix-1;

    % Isolate correct slice & pixel values by adding the best offsets
    firstslice = 1 + Bslc;
    lastslice = size(img_ref,3) + Bslc;
    IMG = SliceCombine(img, ifull+Bpix, ipart);
    MSK = SliceCombine(msk, ifull+Bpix, ipart);
    IMG = IMG(:,:,firstslice:lastslice);
    MSK = MSK(:,:,firstslice:lastslice);

    % Set necessary variables as outputs
    varargout{1} = Bpix;
    varargout{2} = Bslc;
end

% Determine best result
function B = IsolateBestRS(met)
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
IMG = Norm01(IMG);

% Mostly vectorized version of above code. It is less efficient.
% % Sum full slices + part of previous slice + part of next slice
% for y = 1:size(ifull,2)
%     IMG(:,:,y) = sum(img(:,:,ifull(1,y):ifull(2,y)), 3);
% end
% % Add partial slice data (excluding index 0)
% IMG(:,:,2:end) = IMG(:,:,2:end) + img(:,:,ifull(1,2:end)-1).*reshape(ipart(1,2:end), [1,1,size(ipart(1,2:end),2)]);
% IMG = IMG + img(:,:,ifull(2,:)+1).*reshape(ipart(2,:), [1,1,size(ipart,2)]);
end
end

% Calculate mutual information index for given slice set
function I = MII(A, B)
% Inputs
%        A: dataset 1
%        B: dataset 2
% Outputs
%        I: mutual information index

% Preallocate data
I = zeros(1,size(A,3));

% Calculate mutual information
for i = 1:length(I)
    I(i) = mutualInfo(A, B);
end

I = mean(I);
end

% dimlen = size(img_oldres,1); % find image size
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
