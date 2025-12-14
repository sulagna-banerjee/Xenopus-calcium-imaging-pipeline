%% Regional Step 3 — Create central circular masks (fixed pixel count) from L/R MB/FB masks
% and visualize them on the max-fluorescence slice (per replicate)
%
% Works with outputs from:
%   Regional Step 1: Midbrain_selectedRegionMask.mat, Forebrain_selectedRegionMask.mat
%   Regional Step 2: Midbrain_Left_RegionMask.mat, Midbrain_Right_RegionMask.mat,
%                    Forebrain_Left_RegionMask.mat, Forebrain_Right_RegionMask.mat
%
% Expected per replicate folder:
%   - final1.tif
%   - Midbrain_Left_RegionMask.mat      (midbrainLeftMask)
%   - Midbrain_Right_RegionMask.mat     (midbrainRightMask)
%   - Forebrain_Left_RegionMask.mat     (forebrainLeftMask)
%   - Forebrain_Right_RegionMask.mat    (forebrainRightMask)
%
% Outputs (saved to results/ so code-only repo is clean):
%   results/regional_step3_central_points/
%     - <replicateID>_center_overlay.png
%     - <replicateID>_ROI_CentralMasks.mat  (centralMasks, centers, fluoroData, maxSliceIdx)
%
% Folder structure (same as Regional Step 1 & 2):
% project_root/
% ├── data/
% │   ├── control/<replicate_id>/...
% │   └── test/<replicate_id>/...
% └── results/
%     └── regional_step3_central_points/

%% ============================
% 0) USER CONFIGURATION
% ============================

projectRoot = pwd;

dataDir    = fullfile(projectRoot, "data");
controlDir = fullfile(dataDir, "control");
testDir    = fullfile(dataDir, "test");

outputDir = fullfile(projectRoot, "results", "regional_step3_central_points");
if ~exist(outputDir, "dir"), mkdir(outputDir); end

targetPixelCount = 500;  % central mask pixel count for each region

% Provide replicate lists OR auto-discover
controlReplicates = {};      % e.g., {"rep1","rep2"}
testReplicates    = {};      % e.g., {"repA","repB"}
autoDiscoverReplicates = true;

% If true, skip replicate if outputs already exist
skipIfOutputsExist = true;

%% ============================
% 1) RESOLVE REPLICATE LISTS
% ============================

if autoDiscoverReplicates
    if isempty(controlReplicates), controlReplicates = listSubfolders(controlDir); end
    if isempty(testReplicates),    testReplicates    = listSubfolders(testDir);    end
end

allReplicates = [controlReplicates, testReplicates];
allGroups     = [repmat({"control"}, 1, numel(controlReplicates)), repmat({"test"}, 1, numel(testReplicates))];

disp("Regional Step 3 — Central masks + overlay");
disp("Controls: " + numel(controlReplicates) + " | Tests: " + numel(testReplicates));
disp("Output: " + outputDir);

%% ============================
% 2) MAIN LOOP
% ============================

for i = 1:numel(allReplicates)

    replicateID = allReplicates{i};
    groupLabel  = allGroups{i};

    if strcmpi(groupLabel, "control")
        replicateFolder = fullfile(controlDir, replicateID);
    else
        replicateFolder = fullfile(testDir, replicateID);
    end

    disp("--------------------------------------------------");
    disp("Processing replicate: " + replicateID + " (" + groupLabel + ")");
    disp("Folder: " + replicateFolder);

    outPng = fullfile(outputDir, replicateID + "_center_overlay.png");
    outMat = fullfile(outputDir, replicateID + "_ROI_CentralMasks.mat");
    if skipIfOutputsExist && isfile(outPng) && isfile(outMat)
        disp("Outputs already exist — skipping.");
        continue;
    end

    % Inputs
    tifPath = fullfile(replicateFolder, "final1.tif");
    if ~isfile(tifPath)
        warning("Missing final1.tif for %s — skipping.", replicateID);
        continue;
    end

    maskPaths = struct( ...
        "Mid_L",  fullfile(replicateFolder, "Midbrain_Left_RegionMask.mat"), ...
        "Mid_R",  fullfile(replicateFolder, "Midbrain_Right_RegionMask.mat"), ...
        "Fore_L", fullfile(replicateFolder, "Forebrain_Left_RegionMask.mat"), ...
        "Fore_R", fullfile(replicateFolder, "Forebrain_Right_RegionMask.mat") ...
    );

    if ~isfile(maskPaths.Mid_L) || ~isfile(maskPaths.Mid_R) || ~isfile(maskPaths.Fore_L) || ~isfile(maskPaths.Fore_R)
        warning("Missing one or more L/R region masks for %s — run Regional Step 2 first.", replicateID);
        continue;
    end

    % Load max fluorescence slice (memory-safe: does not keep full stack)
    [maxSlice, maxSliceIdx] = loadMaxFluorescenceSlice(tifPath);

    % Load region masks
    midbrainLeftMask   = loadMask(maskPaths.Mid_L,  "midbrainLeftMask");
    midbrainRightMask  = loadMask(maskPaths.Mid_R,  "midbrainRightMask");
    forebrainLeftMask  = loadMask(maskPaths.Fore_L, "forebrainLeftMask");
    forebrainRightMask = loadMask(maskPaths.Fore_R, "forebrainRightMask");

    % Initialize outputs fresh each replicate
    centralMasks = struct();
    centers      = struct();
    fluoroData   = struct();

    % Central circular masks
    [centralMasks.Mid_L,  centers.Mid_L]  = getGeometricCircularMask(midbrainLeftMask,   targetPixelCount);
    [centralMasks.Mid_R,  centers.Mid_R]  = getGeometricCircularMask(midbrainRightMask,  targetPixelCount);
    [centralMasks.Fore_L, centers.Fore_L] = getGeometricCircularMask(forebrainLeftMask,  targetPixelCount);
    [centralMasks.Fore_R, centers.Fore_R] = getGeometricCircularMask(forebrainRightMask, targetPixelCount);

    % Fluorescence extraction on max slice (omitnan-safe)
    fluoroData.Mid_L  = mean(maxSlice(centralMasks.Mid_L),  "omitnan");
    fluoroData.Mid_R  = mean(maxSlice(centralMasks.Mid_R),  "omitnan");
    fluoroData.Fore_L = mean(maxSlice(centralMasks.Fore_L), "omitnan");
    fluoroData.Fore_R = mean(maxSlice(centralMasks.Fore_R), "omitnan");

    % Visualization overlay
    regionMasks = {midbrainLeftMask, midbrainRightMask, forebrainLeftMask, forebrainRightMask};
    fieldNames  = {"Mid_L", "Mid_R", "Fore_L", "Fore_R"};

    colors = {
        [230 159   0]/255;   % Left MB  (orange)
        [ 86 180 233]/255;   % Right MB (sky blue)
        [  0 158 115]/255;   % Left FB  (bluish green)
        [240 228  66]/255    % Right FB (yellow)
    };
    labels = {"Left MB", "Right MB", "Left FB", "Right FB"};

    fig = figure("Visible", "off");
    imshow(maxSlice, []); hold on;

    h = gobjects(1, numel(regionMasks));
    for j = 1:numel(regionMasks)
        contour(regionMasks{j}, [0.5 0.5], "LineWidth", 2, "LineColor", colors{j});

        cm = centralMasks.(fieldNames{j});
        [yC, xC] = find(cm);

        fprintf("  %s central size: %d pixels\n", fieldNames{j}, nnz(cm));

        if isempty(xC) || isempty(yC)
            warning("%s central mask is empty for replicate %s", labels{j}, replicateID);
            continue;
        end

        h(j) = scatter(xC, yC, 10, colors{j}, "filled");
    end

    valid = isgraphics(h);
    if any(valid)
        legend(h(valid), labels(valid), "Location", "eastoutside");
    end
    title("Geometric Centered ROIs (" + targetPixelCount + " px) — " + replicateID);

    saveas(fig, outPng);
    close(fig);

    % Save MAT output
    save(outMat, "centralMasks", "centers", "fluoroData", "maxSliceIdx", "-v7.3");

    fprintf("✅ Saved overlay: %s\n", outPng);
    fprintf("✅ Saved MAT:     %s\n", outMat);
end

disp("Regional Step 3 complete.");

%% ============================================================
% FUNCTIONS
%% ============================================================

function reps = listSubfolders(parentDir)
    reps = {};
    if ~exist(parentDir, "dir")
        warning("Directory not found: %s", parentDir);
        return;
    end
    d = dir(parentDir);
    isGood = [d.isdir] & ~ismember({d.name}, {".",".."});
    reps = sort({d(isGood).name});
end

function mask = loadMask(matPath, varName)
    S = load(matPath, varName);
    assert(isfield(S, varName), "File %s missing variable '%s'.", matPath, varName);
    mask = logical(S.(varName));
end

function [maxSlice, maxIdx] = loadMaxFluorescenceSlice(tifPath)
    info = imfinfo(tifPath);
    nSlices = numel(info);

    meanFluoro = zeros(nSlices, 1);
    for k = 1:nSlices
        img = im2double(imread(tifPath, k));
        meanFluoro(k) = mean(img, "all");
    end

    [~, maxIdx] = max(meanFluoro);
    maxSlice = mat2gray(im2double(imread(tifPath, maxIdx)));
end

function [circleMask, center] = getGeometricCircularMask(regionMask, targetPixelCount)
    regionMask = bwareafilt(logical(regionMask), 1); % largest connected component
    stats = regionprops(regionMask, "Centroid");
    assert(~isempty(stats), "Region mask is empty after filtering.");
    center = stats.Centroid; % [x, y]

    r = sqrt(targetPixelCount / pi);
    [X, Y] = meshgrid(1:size(regionMask,2), 1:size(regionMask,1));
    circle = ((X - center(1)).^2 + (Y - center(2)).^2) <= r^2;

    circleMask = circle & regionMask;

    % Trim if needed (keep densest interior pixels)
    if nnz(circleMask) > targetPixelCount
        d = bwdist(~circleMask);
        [~, idx] = maxk(d(:), targetPixelCount);
        trimmed = false(size(regionMask));
        trimmed(idx) = true;
        circleMask = trimmed;
    end
end
