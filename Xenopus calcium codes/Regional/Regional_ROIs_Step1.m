%% Regional Step 1 — Draw midbrain/forebrain masks and extract regional ROI label maps
%
% Public/GitHub-ready script (code-only release).
%
% This script assumes you already ran the whole-brain ROI step (Whole-brain Step 1),
% which produces, per replicate folder:
%   - WholeBrain_selectedRegionMask.mat            (variable: wholeBrainMask)
%   - WholeBrain_Voronoi_Refined_ROIs.mat          (variable: granules_voronoi_wholeBrain)
%
% It then:
%   1) Loads final1.tif for that replicate
%   2) Finds the max-fluorescence slice
%   3) Overlays whole-brain perimeter for guidance
%   4) Prompts user to draw Midbrain and Forebrain regions (unless masks already exist)
%   5) Applies masks to whole-brain ROI label map to generate regional ROI label maps
%
% -----------------------------
% Suggested folder structure:
%
% project_root/
% ├── scripts/
% │   └── regional_step1_draw_masks_and_extract_rois.m
% ├── data/
% │   ├── control/
% │   │   └── <replicate_id>/
% │   │       ├── final1.tif
% │   │       ├── WholeBrain_selectedRegionMask.mat
% │   │       └── WholeBrain_Voronoi_Refined_ROIs.mat
% │   └── test/
% │       └── <replicate_id>/
% │           ├── final1.tif
% │           ├── WholeBrain_selectedRegionMask.mat
% │           └── WholeBrain_Voronoi_Refined_ROIs.mat
%
% Outputs written into each replicate folder:
%   - Midbrain_selectedRegionMask.mat              (midbrainMask)
%   - Forebrain_selectedRegionMask.mat             (forebrainMask)
%   - Midbrain_Voronoi_Refined_ROIs.mat            (granules_voronoi_midbrain)
%   - Forebrain_Voronoi_Refined_ROIs.mat           (granules_voronoi_forebrain)
%
% -----------------------------
% Requirements:
%   - Image Processing Toolbox (imfreehand/createMask/bwperim/imoverlay)
%   - Parallel Computing Toolbox (optional; used for TIFF loading with parfor)
%
% Author: <Author/Lab>
% Associated manuscript: <Paper title/DOI>

%% ============================
% 0) USER CONFIGURATION
% ============================

projectRoot = pwd;

dataDir    = fullfile(projectRoot, "data");
controlDir = fullfile(dataDir, "control");
testDir    = fullfile(dataDir, "test");

% Provide replicate lists OR enable auto-discovery
controlReplicates = {};    % e.g., {"rep1","rep2"}
testReplicates    = {};    % e.g., {"repA","repB"}

autoDiscoverReplicates = true;

%% ============================
% 1) PARALLEL POOL CLEANUP (OPTIONAL)
% ============================

if ~isempty(gcp("nocreate"))
    disp("Closing existing parallel pool to free memory...");
    delete(gcp("nocreate"));
end

%% ============================
% 2) RESOLVE REPLICATE LISTS
% ============================

if autoDiscoverReplicates
    if isempty(controlReplicates)
        controlReplicates = listSubfolders(controlDir);
        disp("Auto-discovered control replicates: " + numel(controlReplicates));
    end
    if isempty(testReplicates)
        testReplicates = listSubfolders(testDir);
        disp("Auto-discovered test replicates: " + numel(testReplicates));
    end
end

allReplicates = [controlReplicates, testReplicates];
allGroups     = [repmat({"control"}, 1, numel(controlReplicates)), repmat({"test"}, 1, numel(testReplicates))];

%% ============================
% 3) MAIN LOOP
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

    % Required inputs from Whole-brain Step 1
    wholeBrainMaskFile = fullfile(replicateFolder, "WholeBrain_selectedRegionMask.mat");
    wholeBrainROIsFile = fullfile(replicateFolder, "WholeBrain_Voronoi_Refined_ROIs.mat");
    tifFilePath        = fullfile(replicateFolder, "final1.tif");

    assert(exist(wholeBrainMaskFile, "file") == 2, "Missing: %s", wholeBrainMaskFile);
    assert(exist(wholeBrainROIsFile, "file") == 2, "Missing: %s", wholeBrainROIsFile);
    assert(exist(tifFilePath, "file") == 2,        "Missing: %s", tifFilePath);

    % Regional mask files
    midbrainMaskFile  = fullfile(replicateFolder, "Midbrain_selectedRegionMask.mat");
    forebrainMaskFile = fullfile(replicateFolder, "Forebrain_selectedRegionMask.mat");

    % Load whole-brain mask + ROI labels
    disp("Loading whole-brain mask + ROI labels...");
    S1 = load(wholeBrainMaskFile, "wholeBrainMask");
    S2 = load(wholeBrainROIsFile, "granules_voronoi_wholeBrain");

    wholeBrainMask = S1.wholeBrainMask;
    granules_voronoi_wholeBrain = S2.granules_voronoi_wholeBrain;

    % Ensure ROI label map is 3D for slice-wise multiplication
    if ndims(granules_voronoi_wholeBrain) == 2
        granules_voronoi_wholeBrain = repmat(granules_voronoi_wholeBrain, [1 1 1]);
    end

    % If masks exist, load them; otherwise ask user to draw
    if exist(midbrainMaskFile, "file") == 2 && exist(forebrainMaskFile, "file") == 2
        disp("Regional masks already exist — loading and skipping drawing.");
        load(midbrainMaskFile, "midbrainMask");
        load(forebrainMaskFile, "forebrainMask");
    else
        % Load TIFF and find max fluorescence slice for drawing guidance
        disp("Loading TIFF stack for drawing...");
        info = imfinfo(tifFilePath);
        nSlices = numel(info);

        imageStack = zeros(info(1).Height, info(1).Width, nSlices, "double");

        % Parallel load if possible; fall back to for-loop if parfor unavailable
        try
            parfor k = 1:nSlices
                imageStack(:,:,k) = im2double(imread(tifFilePath, k));
            end
        catch
            for k = 1:nSlices
                imageStack(:,:,k) = im2double(imread(tifFilePath, k));
            end
        end

        % Find max fluorescence slice
        disp("Finding max-fluorescence slice...");
        maxIdx = 1;
        maxVal = -inf;
        for s = 1:nSlices
            m = mean(imageStack(:,:,s), "all");
            if m > maxVal
                maxVal = m;
                maxIdx = s;
            end
        end
        disp("Max fluorescence slice index: " + maxIdx);

        maxSlice = imageStack(:,:,maxIdx);

        % Overlay whole-brain perimeter for guidance
        wholeBrainPerimeter = bwperim(wholeBrainMask);
        maxSliceOverlay = imoverlay(mat2gray(maxSlice), wholeBrainPerimeter, [1 0 0]);

        % Draw midbrain
        figure;
        imshow(maxSliceOverlay, []);
        title("Draw MIDBRAIN region — " + replicateID);
        hMid = imfreehand;
        midbrainMask = createMask(hMid) > 0;
        save(midbrainMaskFile, "midbrainMask");
        close;

        % Draw forebrain
        figure;
        imshow(maxSliceOverlay, []);
        title("Draw FOREBRAIN region — " + replicateID);
        hFore = imfreehand;
        forebrainMask = createMask(hFore) > 0;
        save(forebrainMaskFile, "forebrainMask");
        close;

        clear imageStack;
        disp("Regional masks drawn and saved.");
    end

    % Output ROI label maps
    midbrainROIsFile  = fullfile(replicateFolder, "Midbrain_Voronoi_Refined_ROIs.mat");
    forebrainROIsFile = fullfile(replicateFolder, "Forebrain_Voronoi_Refined_ROIs.mat");

    if exist(midbrainROIsFile, "file") == 2 && exist(forebrainROIsFile, "file") == 2
        disp("Regional ROI label maps already exist — skipping ROI extraction.");
        continue;
    end

    % Apply masks slice-wise to ROI labels
    disp("Applying midbrain mask to whole-brain ROI label map...");
    granules_voronoi_midbrain = granules_voronoi_wholeBrain;
    for s = 1:size(granules_voronoi_wholeBrain, 3)
        granules_voronoi_midbrain(:,:,s) = granules_voronoi_wholeBrain(:,:,s) .* midbrainMask;
    end
    save(midbrainROIsFile, "granules_voronoi_midbrain", "-v7.3");

    disp("Applying forebrain mask to whole-brain ROI label map...");
    granules_voronoi_forebrain = granules_voronoi_wholeBrain;
    for s = 1:size(granules_voronoi_wholeBrain, 3)
        granules_voronoi_forebrain(:,:,s) = granules_voronoi_wholeBrain(:,:,s) .* forebrainMask;
    end
    save(forebrainROIsFile, "granules_voronoi_forebrain", "-v7.3");

    % Cleanup
    clear wholeBrainMask midbrainMask forebrainMask granules_voronoi_wholeBrain;
    disp("Saved regional ROI label maps.");
end

%% ============================
% 4) CLOSE PARALLEL POOL (OPTIONAL)
% ============================

if ~isempty(gcp("nocreate"))
    disp("Closing parallel pool...");
    delete(gcp("nocreate"));
end

disp("Regional segmentation completed.");

%% ============================================================
% Helper: list replicate folders in a directory
%% ============================================================
function reps = listSubfolders(parentDir)
    reps = {};
    if ~exist(parentDir, "dir")
        warning("Directory not found: %s", parentDir);
        return;
    end
    d = dir(parentDir);
    isGood = [d.isdir] & ~ismember({d.name}, {".",".."});
    reps = {d(isGood).name};
    reps = sort(reps);
end
