%% Step 2 — Collect raw ROI fluorescence traces for whole-brain comparisons
%
% Public/GitHub-ready script (code-only release).
% Users supply their own replicate folder names and data.
%
% -----------------------------
% Expected folder structure:
%
% project_root/
% ├── scripts/
% │   └── step2_collect_raw_fluorescence_whole_brain.m
% ├── data/
% │   ├── control/
% │   │   └── <replicate_id>/Filtered_Fluorescence_Data.mat
% │   └── test/
% │       └── <replicate_id>/Filtered_Fluorescence_Data.mat
% └── results/
%
% Each Filtered_Fluorescence_Data.mat must contain:
%   fluorescenceData   [nROIs x nFrames]
%
% -----------------------------
% Outputs:
%   results/step2_raw_trace_collection/
%     - Control_Group_Data.mat
%     - Test_Group_Data.mat
%     - Control_Summary.csv
%     - Test_Summary.csv
%
% Notes:
%   This script *collects raw traces only* (no peak detection here).
%   Whole-brain relevance depends on how Filtered_Fluorescence_Data.mat was created
%   (e.g., ROI traces derived from a whole-brain mask in Step 1).
%
% Author: <Author/Lab>
% Associated manuscript: <Paper title/DOI>

%% ============================
% 0) USER CONFIGURATION
% ============================

projectRoot = pwd;

% GitHub-clean default locations
dataDir       = fullfile(projectRoot, "data");
controlDir    = fullfile(dataDir, "control");
testDir       = fullfile(dataDir, "test");

outputDir     = fullfile(projectRoot, "results", "step2_raw_trace_collection");
if ~exist(outputDir, "dir")
    mkdir(outputDir);
end

% Users: provide replicate folder names here (or leave empty and use auto-discovery below)
controlReplicates = {};   % e.g., {"rep1","rep2","rep3"}
testReplicates    = {};   % e.g., {"repA","repB"}

% Optional: auto-discover replicate folders if lists are left empty
autoDiscoverReplicates = true;

%% ============================
% 1) RESOLVE REPLICATE LISTS
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

% Basic checks
if isempty(controlReplicates)
    warning("No control replicates provided/found. Control outputs will be empty.");
end
if isempty(testReplicates)
    warning("No test replicates provided/found. Test outputs will be empty.");
end

%% ============================
% 2) INITIALISE STRUCTURES
% ============================

% Preallocate to max possible size; we'll trim later if some replicates fail/skip
controlData = repmat(struct( ...
    "fluorescenceData", [], ...
    "replicateID", "", ...
    "groupLabel", "", ...
    "roiCount", 0, ...
    "frameCount", 0), 1, numel(controlReplicates));

testData = repmat(struct( ...
    "fluorescenceData", [], ...
    "replicateID", "", ...
    "groupLabel", "", ...
    "roiCount", 0, ...
    "frameCount", 0), 1, numel(testReplicates));

controlSummary = {};
testSummary    = {};

controlKeep = false(1, numel(controlReplicates));
testKeep    = false(1, numel(testReplicates));

%% ============================
% 3) LOAD CONTROL DATA
% ============================

disp("Loading control group data...");
for i = 1:numel(controlReplicates)
    repID = controlReplicates{i};

    try
        fluorescenceData = loadRawFluorescence(controlDir, repID);

        if isempty(fluorescenceData) || ndims(fluorescenceData) ~= 2
            warning("Skipping control replicate %s: invalid fluorescenceData (empty or not 2D).", repID);
            continue;
        end

        controlData(i).fluorescenceData = fluorescenceData;
        controlData(i).replicateID      = repID;
        controlData(i).groupLabel       = "control";
        controlData(i).roiCount         = size(fluorescenceData, 1);
        controlData(i).frameCount       = size(fluorescenceData, 2);

        controlSummary(end+1, :) = {repID, "control", size(fluorescenceData, 1), size(fluorescenceData, 2)};
        controlKeep(i) = true;

    catch ME
        warning("Failed to load control replicate %s: %s", repID, ME.message);
    end
end

%% ============================
% 4) LOAD TEST DATA
% ============================

disp("Loading test group data...");
for i = 1:numel(testReplicates)
    repID = testReplicates{i};

    try
        fluorescenceData = loadRawFluorescence(testDir, repID);

        if isempty(fluorescenceData) || ndims(fluorescenceData) ~= 2
            warning("Skipping test replicate %s: invalid fluorescenceData (empty or not 2D).", repID);
            continue;
        end

        testData(i).fluorescenceData = fluorescenceData;
        testData(i).replicateID      = repID;
        testData(i).groupLabel       = "test";
        testData(i).roiCount         = size(fluorescenceData, 1);
        testData(i).frameCount       = size(fluorescenceData, 2);

        testSummary(end+1, :) = {repID, "test", size(fluorescenceData, 1), size(fluorescenceData, 2)};
        testKeep(i) = true;

    catch ME
        warning("Failed to load test replicate %s: %s", repID, ME.message);
    end
end

% Trim out skipped entries (keeps files cleaner)
controlData = controlData(controlKeep);
testData    = testData(testKeep);

%% ============================
% 5) SAVE OUTPUTS
% ============================

disp("Saving group-based raw fluorescence data...");
save(fullfile(outputDir, "Control_Group_Data.mat"), "controlData", "-v7.3");
save(fullfile(outputDir, "Test_Group_Data.mat"), "testData", "-v7.3");
disp("Raw data saved successfully.");

% Save summary CSVs
summaryHeader = {"ReplicateID", "Group", "ROICount", "FrameCount"};

controlSummaryTable = cell2table(controlSummary, "VariableNames", summaryHeader);
testSummaryTable    = cell2table(testSummary, "VariableNames", summaryHeader);

writetable(controlSummaryTable, fullfile(outputDir, "Control_Summary.csv"));
writetable(testSummaryTable, fullfile(outputDir, "Test_Summary.csv"));

disp("Summary CSV files saved.");
disp("Step 2 complete.");

%% ============================================================
% Helper: Load raw fluorescence data from one replicate
%% ============================================================
function fluorescenceData = loadRawFluorescence(groupDir, replicateID)
    disp("Loading replicate: " + replicateID);

    filePath = fullfile(groupDir, replicateID, "Filtered_Fluorescence_Data.mat");
    disp("Loading file: " + filePath);

    if ~exist(filePath, "file")
        error("Missing file: %s", filePath);
    end

    S = load(filePath, "fluorescenceData");
    if ~isfield(S, "fluorescenceData")
        error("File does not contain variable 'fluorescenceData': %s", filePath);
    end

    fluorescenceData = S.fluorescenceData;
end

%% ============================================================
% Helper: List subfolders (replicate IDs) within a directory
%% ============================================================
function reps = listSubfolders(parentDir)
    reps = {};
    if ~exist(parentDir, "dir")
        warning("Directory not found: %s", parentDir);
        return;
    end

    d = dir(parentDir);
    isGood = [d.isdir] & ~ismember({d.name}, {".", ".."});
    reps = {d(isGood).name};

    % Sort for reproducibility
    reps = sort(reps);
end
