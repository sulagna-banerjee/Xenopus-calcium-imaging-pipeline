%% Regional Step 2 — Manually split midbrain and forebrain into left/right masks
%
% Public/GitHub-ready script (code-only release).
%
% Prerequisites (per replicate folder):
%   - final1.tif
%   - Midbrain_selectedRegionMask.mat   (midbrainMask)
%   - Forebrain_selectedRegionMask.mat  (forebrainMask)
%
% Outputs (per replicate folder):
%   - Midbrain_Left_RegionMask.mat      (midbrainLeftMask)
%   - Midbrain_Right_RegionMask.mat     (midbrainRightMask)
%   - Forebrain_Left_RegionMask.mat     (forebrainLeftMask)
%   - Forebrain_Right_RegionMask.mat    (forebrainRightMask)
%
% Suggested folder structure:
% project_root/
% ├── scripts/
% │   └── regional_step2_split_mb_fb_left_right.m
% └── data/
%     ├── control/<replicate_id>/...
%     └── test/<replicate_id>/...
%
% Notes:
%   - Uses max-fluorescence slice for drawing
%   - Overlays MB/FB contours for guidance
%   - Requires Image Processing Toolbox (imfreehand/createMask/contour)
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
controlReplicates = {};  % e.g., {"rep1","rep2"}
testReplicates    = {};  % e.g., {"repA","repB"}
autoDiscoverReplicates = true;

% If true, skip drawing when the 4 L/R mask files already exist
skipIfAlreadyExists = true;

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

    % Required inputs
    tifFilePath = fullfile(replicateFolder, "final1.tif");
    midbrainMaskFile  = fullfile(replicateFolder, "Midbrain_selectedRegionMask.mat");
    forebrainMaskFile = fullfile(replicateFolder, "Forebrain_selectedRegionMask.mat");

    assert(isfile(tifFilePath),        "Missing TIFF file: %s", tifFilePath);
    assert(isfile(midbrainMaskFile),   "Missing midbrain mask: %s", midbrainMaskFile);
    assert(isfile(forebrainMaskFile),  "Missing forebrain mask: %s", forebrainMaskFile);

    % Outputs
    outMBL = fullfile(replicateFolder, "Midbrain_Left_RegionMask.mat");
    outMBR = fullfile(replicateFolder, "Midbrain_Right_RegionMask.mat");
    outFBL = fullfile(replicateFolder, "Forebrain_Left_RegionMask.mat");
    outFBR = fullfile(replicateFolder, "Forebrain_Right_RegionMask.mat");

    if skipIfAlreadyExists && isfile(outMBL) && isfile(outMBR) && isfile(outFBL) && isfile(outFBR)
        disp("Left/right masks already exist — skipping drawing.");
        continue;
    end

    % Load TIFF stack and find max fluorescence slice
    disp("Loading TIFF stack...");
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

    maxSlice = mat2gray(imageStack(:,:,maxIdx));

    % Load existing MB/FB masks
    Smb = load(midbrainMaskFile, "midbrainMask");
    Sfb = load(forebrainMaskFile, "forebrainMask");
    midbrainMask  = Smb.midbrainMask > 0;
    forebrainMask = Sfb.forebrainMask > 0;

    %% ---- Draw LEFT/RIGHT Midbrain ----
    fig1 = figure;
    imshow(maxSlice, []);
    hold on;
    contour(midbrainMask, [0.5 0.5], "r", "LineWidth", 2);

    title("Draw LEFT Midbrain — " + replicateID);
    hMidLeft = imfreehand;
    midbrainLeftMask = createMask(hMidLeft) > 0;
    save(outMBL, "midbrainLeftMask");
    disp("Saved left midbrain mask.");

    title("Draw RIGHT Midbrain — " + replicateID);
    hMidRight = imfreehand;
    midbrainRightMask = createMask(hMidRight) > 0;
    save(outMBR, "midbrainRightMask");
    disp("Saved right midbrain mask.");

    close(fig1);

    %% ---- Draw LEFT/RIGHT Forebrain ----
    fig2 = figure;
    imshow(maxSlice, []);
    hold on;
    contour(forebrainMask, [0.5 0.5], "b", "LineWidth", 2);

    title("Draw LEFT Forebrain — " + replicateID);
    hForeLeft = imfreehand;
    forebrainLeftMask = createMask(hForeLeft) > 0;
    save(outFBL, "forebrainLeftMask");
    disp("Saved left forebrain mask.");

    title("Draw RIGHT Forebrain — " + replicateID);
    hForeRight = imfreehand;
    forebrainRightMask = createMask(hForeRight) > 0;
    save(outFBR, "forebrainRightMask");
    disp("Saved right forebrain mask.");

    close(fig2);

    clear imageStack midbrainMask forebrainMask;
end

disp("Splitting of midbrain and forebrain into left/right regions completed.");

%% ============================
% 4) CLOSE PARALLEL POOL (OPTIONAL)
% ============================

if ~isempty(gcp("nocreate"))
    disp("Closing parallel pool...");
    delete(gcp("nocreate"));
end

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
