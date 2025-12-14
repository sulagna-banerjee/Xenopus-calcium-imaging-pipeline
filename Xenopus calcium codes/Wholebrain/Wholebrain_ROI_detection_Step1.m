%% Step 1 — Whole-brain ROI detection (single replicate)
% Purpose:
%   Load a TIFF time-series stack (final1.tif), draw a whole-brain mask,
%   apply the mask across the stack, run CalciSeg ROI detection, and save
%   outputs + a visualization overlay for QC.
%
% Inputs (expected):
%   - A replicate subfolder containing: final1.tif
%   - CalciSeg on the MATLAB path (or set CALCISEG_DIR below)
%
% Outputs (written to replicateOutputFolder):
%   - Loaded_ImageStack.mat
%   - Max_Fluorescence_Slice.mat
%   - WholeBrain_selectedRegionMask.mat
%   - Masked_WholeBrain_ImageStack.mat
%   - WholeBrain_Voronoi_Refined_ROIs.mat
%   - WholeBrain_ROIs_Voronoi_Refined_Visualization_Overlay.png
%
% Notes:
%   - This script uses interactive ROI drawing (imfreehand).
%   - TIFF loading is parallelized using parfor (Parallel Computing Toolbox).
%
% Suggested repo structure:
%   repo/
%     scripts/
%     data/               (optional example data; don’t upload raw large files)
%     results/            (or user-defined output)
%     external/CalciSeg/  (optional submodule / vendor drop)
%
% Author: (Your name / lab)
% Associated manuscript: (Paper title / DOI placeholder)

%% ============================
% 0) User-configurable settings
% ============================

% ---- (A) CalciSeg location ----
% Option 1: keep CalciSeg as a subfolder in the repo (recommended)
%   e.g., external/CalciSeg
% Option 2: point to a local installation (set via environment variable / config)
CALCISEG_DIR = fullfile(pwd, "external", "CalciSeg"); % <-- recommended default

% If you don't keep CalciSeg in-repo, set CALCISEG_DIR to your local folder,
% or leave empty and ensure CalciSeg is already on the MATLAB path.
% CALCISEG_DIR = "PATH_TO_CALCISEG";

% ---- (B) Replicate selection ----
replicateSubfolder   = "20250219i";          % replicate folder name
replicateOutputName  = replicateSubfolder;   % output folder name (usually same)

% ---- (C) Input/output base directories ----
% Use relative paths for GitHub friendliness.
baseFolder   = fullfile(pwd, "data");     % contains replicateSubfolder/final1.tif
outputFolder = fullfile(pwd, "results");  % outputs will be placed here

% ---- (D) File name ----
tifFilename = "final1.tif";

%% ============================
% 1) Safety + path setup
% ============================

% Close parallel pool if it exists (frees memory)
poolObj = gcp("nocreate");
if ~isempty(poolObj)
    disp("Closing existing parallel pool to free memory...");
    delete(poolObj);
end

% Add CalciSeg to path (if provided)
if ~isempty(CALCISEG_DIR)
    if isfolder(CALCISEG_DIR)
        addpath(genpath(CALCISEG_DIR));
        disp("Added CalciSeg to MATLAB path:");
        disp("  " + string(CALCISEG_DIR));
    else
        error("CALCISEG_DIR does not exist: %s", CALCISEG_DIR);
    end
end

% Confirm CalciSeg is found
whichCalciSeg = which("CalciSeg");
if isempty(whichCalciSeg)
    error("CalciSeg not found on MATLAB path. Please set CALCISEG_DIR correctly.");
else
    disp("CalciSeg found at:");
    disp("  " + string(whichCalciSeg));
end

%% ============================
% 2) Resolve paths + output dirs
% ============================

stackFolder  = fullfile(baseFolder, replicateSubfolder);
tifFilePath  = fullfile(stackFolder, tifFilename);

replicateOutputFolder = fullfile(outputFolder, replicateOutputName);
if ~exist(replicateOutputFolder, "dir")
    mkdir(replicateOutputFolder);
    disp("Created output folder for replicate:");
    disp("  " + string(replicateOutputFolder));
end

%% ============================
% 3) Load TIFF stack (parallel)
% ============================

disp("Loading TIFF file for current replicate...");

if ~exist(tifFilePath, "file")
    error("TIFF file not found: %s", tifFilePath);
end

info      = imfinfo(tifFilePath);
numImages = numel(info);

% Preallocate stack
imageStack = zeros(info(1).Height, info(1).Width, numImages, "double");

% Parallel read for speed
parfor k = 1:numImages
    imageStack(:,:,k) = im2double(imread(tifFilePath, k));
end

disp("Image stack loaded successfully.");

% Save loaded stack (optional but useful for checkpointing)
save(fullfile(replicateOutputFolder, "Loaded_ImageStack.mat"), "imageStack", "-v7.3");
disp("Saved Loaded_ImageStack.mat");

%% ============================
% 4) Find max-fluorescence slice
% ============================

disp("Finding slice with maximum mean fluorescence...");
maxFluorescence = -inf;
maxFluorescenceSliceIndex = 1;

for sliceIndex = 1:numImages
    meanFluorescence = mean(imageStack(:,:,sliceIndex), "all");
    if meanFluorescence > maxFluorescence
        maxFluorescence = meanFluorescence;
        maxFluorescenceSliceIndex = sliceIndex;
    end
end

disp("Max fluorescence slice index: " + string(maxFluorescenceSliceIndex));

save(fullfile(replicateOutputFolder, "Max_Fluorescence_Slice.mat"), ...
     "maxFluorescenceSliceIndex");
disp("Saved Max_Fluorescence_Slice.mat");

maxFluorescenceSlice = imageStack(:,:,maxFluorescenceSliceIndex);

%% ============================
% 5) Draw whole-brain mask
% ============================

disp("Draw the whole brain region on the max-fluorescence slice...");

figure;
imshow(maxFluorescenceSlice, []);
title("Draw Whole Brain Region");

hWholeBrain = imfreehand;                     % interactive freehand ROI
wholeBrainMask = createMask(hWholeBrain) > 0; % ensure binary

save(fullfile(replicateOutputFolder, "WholeBrain_selectedRegionMask.mat"), ...
     "wholeBrainMask");
disp("Saved WholeBrain_selectedRegionMask.mat");

%% ============================
% 6) Apply mask to full stack
% ============================

disp("Applying whole-brain mask to the full stack...");

regionMask3D_wholeBrain = repmat(wholeBrainMask, [1, 1, size(imageStack, 3)]);
maskedImageStack_wholeBrain = imageStack .* regionMask3D_wholeBrain;

save(fullfile(replicateOutputFolder, "Masked_WholeBrain_ImageStack.mat"), ...
     "maskedImageStack_wholeBrain", "-v7.3");
disp("Saved Masked_WholeBrain_ImageStack.mat");

% Free memory
clear imageStack wholeBrainMask;
disp("Cleared original imageStack and wholeBrainMask from memory.");

%% ============================
% 7) ROI detection (CalciSeg)
% ============================

disp("Starting ROI detection for whole brain...");

try
    [granules_voronoi_wholeBrain, ~] = CalciSeg(maskedImageStack_wholeBrain, ...
        "projection_method", "corr", ...       % temporal correlation method
        "init_seg_method", "voronoi", ...
        "n_rep", 1, ...
        "limitPixCount", [10, 2000], ...       % allow small + large signals
        "regmax_method", "filtered", ...       % smoothed intensity for seed detection
        "refinement_method", "corr");          % refine based on temporal correlation

    disp("ROI detection completed.");

    % Expand ROI labels to 3D if CalciSeg returned a 2D label image
    if ndims(granules_voronoi_wholeBrain) == 2
        granules_voronoi_wholeBrain = repmat(granules_voronoi_wholeBrain, ...
            [1, 1, size(regionMask3D_wholeBrain, 3)]);
        disp("Expanded granules_voronoi_wholeBrain to 3D.");
    end

    % Remove ROIs outside the mask
    granules_voronoi_wholeBrain(~regionMask3D_wholeBrain) = 0;
    disp("Removed ROIs outside the whole-brain mask.");

    save(fullfile(replicateOutputFolder, "WholeBrain_Voronoi_Refined_ROIs.mat"), ...
         "granules_voronoi_wholeBrain", "-v7.3");
    disp("Saved WholeBrain_Voronoi_Refined_ROIs.mat");

catch ME
    warning("Error during ROI detection:\n%s", getReport(ME, "extended"));
end

%% ============================
% 8) Visualization overlay (QC)
% ============================

disp("Visualizing detected whole-brain ROIs on max-fluorescence slice...");

maxFluorescenceSlice = maskedImageStack_wholeBrain(:,:,maxFluorescenceSliceIndex);

figure;
imshow(maxFluorescenceSlice, []);
hold on;

uniqueROIs = unique(granules_voronoi_wholeBrain);
uniqueROIs(uniqueROIs == 0) = [];

for roiIdx = 1:numel(uniqueROIs)
    roiLabel = uniqueROIs(roiIdx);
    roiMask  = (granules_voronoi_wholeBrain == roiLabel);

    if size(roiMask, 3) >= maxFluorescenceSliceIndex
        [y, x] = find(roiMask(:,:,maxFluorescenceSliceIndex));

        if ~isempty(y) && ~isempty(x)
            rectangle("Position", [min(x), min(y), max(x)-min(x)+1, max(y)-min(y)+1], ...
                "EdgeColor", "g", "LineWidth", 1);
            text(mean(x), mean(y), num2str(roiLabel), ...
                "Color", "w", "FontSize", 10, "HorizontalAlignment", "center");
        else
            warning("ROI %d has no pixels in the max-fluorescence slice.", roiLabel);
        end
    else
        warning("roiMask has insufficient slices; skipping ROI %d.", roiLabel);
    end
end

title("Detected ROIs in Whole Brain (Green)");
hold off;

saveas(gcf, fullfile(replicateOutputFolder, ...
    "WholeBrain_ROIs_Voronoi_Refined_Visualization_Overlay.png"));
disp("Saved ROI overlay visualization PNG.");

%% ============================
% 9) Cleanup
% ============================

poolObj = gcp("nocreate");
if ~isempty(poolObj)
    disp("Closing parallel pool...");
    delete(poolObj);
end

disp("Processing completed for replicate: " + string(replicateSubfolder));
