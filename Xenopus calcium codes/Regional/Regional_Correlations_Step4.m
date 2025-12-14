%% Regional Step 4 — Group-level correlation between central ROI fluorescence (MB/FB L/R)
%
% GitHub-clean script aligned with Regional Steps 1–3.
%
% This script uses:
%   (A) Regional Step 3 outputs (saved in results/regional_step3_central_points):
%       - <replicateID>_ROI_CentralMasks.mat   (centralMasks)
%
%   (B) Per-replicate ROI label maps (from Regional Step 1, plus optional L/R ROI maps if you have them):
%       Option 1 (recommended, simplest): use whole-brain ROI label map + centralMasks
%           - WholeBrain_Voronoi_Refined_ROIs.mat   (granules_voronoi_wholeBrain)
%       Option 2 (if you have separate L/R ROI maps): specify filenames below.
%
%   (C) Per-replicate fluorescence time series produced earlier in your pipeline:
%       - Filtered_Fluorescence_Data.mat containing:
%           fluorescenceData  (nROIs x nFrames)
%           uniqueROIs        (nROIs x 1) ROI label IDs matching the ROI label map
%
% Outputs:
%   results/regional_step4_correlations/
%     - Pairwise_CentralROI_Correlations.csv
%     - Control_Group_CentralROI_CorrelationMatrix.csv
%     - Test_Group_CentralROI_CorrelationMatrix.csv
%     - Control_Group_CentralROI_CorrelationMatrix.png
%     - Test_Group_CentralROI_CorrelationMatrix.png
%     - <Region1>_vs_<Region2>_Correlation.csv   (one per pair)
%
% Folder structure (same as earlier GitHub-clean steps):
% project_root/
% ├── data/
% │   ├── control/<replicate_id>/...
% │   └── test/<replicate_id>/...
% └── results/
%     ├── regional_step3_central_points/...
%     └── regional_step4_correlations/...

%% ============================
% 0) USER CONFIGURATION
% ============================

projectRoot = pwd;

dataDir    = fullfile(projectRoot, "data");
controlDir = fullfile(dataDir, "control");
testDir    = fullfile(dataDir, "test");

inputFolder9a = fullfile(projectRoot, "results", "regional_step3_central_points");
outputFolder  = fullfile(projectRoot, "results", "regional_step4_correlations");
if ~exist(outputFolder, "dir"), mkdir(outputFolder); end

% Provide replicate lists OR auto-discover
controlReplicates = {};
testReplicates    = {};
autoDiscoverReplicates = true;

% Regions (fieldnames must match Step 3 centralMasks)
regions    = {"Mid_L","Mid_R","Fore_L","Fore_R"};
labelNames = {"Left Midbrain","Right Midbrain","Left Forebrain","Right Forebrain"};

% ---- ROI label map source ----
% Use whole-brain ROI label map and simply sample ROI labels under each central mask.
% This keeps the pipeline consistent and avoids requiring extra L/R ROI-map files.
roiLabelMapFilename = "WholeBrain_Voronoi_Refined_ROIs.mat";
roiLabelMapVarname  = "granules_voronoi_wholeBrain";

% ---- Fluorescence data file ----
fluoroMatFilename = "Filtered_Fluorescence_Data.mat";
fluoroVarName     = "fluorescenceData";
roiIdVarName      = "uniqueROIs";

% Correlation type
corrType = "Pearson";

%% ============================
% 1) RESOLVE REPLICATE LISTS
% ============================

if autoDiscoverReplicates
    if isempty(controlReplicates), controlReplicates = listSubfolders(controlDir); end
    if isempty(testReplicates),    testReplicates    = listSubfolders(testDir);    end
end

allGroups = {"Control","Test"};

% Storage
controlMatrixStack = [];
testMatrixStack    = [];
pairwiseCorrelationData = {};

%% ============================
% 2) PROCESS GROUPS
% ============================

for g = 1:2
    groupName = allGroups{g};
    if strcmpi(groupName, "Control")
        replicates = controlReplicates;
        groupDir   = controlDir;
    else
        replicates = testReplicates;
        groupDir   = testDir;
    end

    for r = 1:numel(replicates)
        replicateID = replicates{r};
        replicateFolder = fullfile(groupDir, replicateID);

        try
            % Step 3 output (central masks)
            roiFile = fullfile(inputFolder9a, replicateID + "_ROI_CentralMasks.mat");
            if ~isfile(roiFile)
                warning("Missing Step 3 file for %s — skipping.", replicateID);
                continue;
            end
            S = load(roiFile, "centralMasks");
            centralMasks = S.centralMasks;

            % Fluorescence time series
            fluoroFile = fullfile(replicateFolder, fluoroMatFilename);
            if ~isfile(fluoroFile)
                warning("Missing fluorescence file for %s — skipping.", replicateID);
                continue;
            end
            F = load(fluoroFile, fluoroVarName, roiIdVarName);
            fluorescenceData = F.(fluoroVarName);   % nROIs x nFrames
            uniqueROIs       = F.(roiIdVarName);    % nROIs x 1 label IDs

            % ROI label map (whole brain)
            roiMapFile = fullfile(replicateFolder, roiLabelMapFilename);
            if ~isfile(roiMapFile)
                warning("Missing ROI label map for %s — skipping.", replicateID);
                continue;
            end
            M = load(roiMapFile, roiLabelMapVarname);
            roiMap = M.(roiLabelMapVarname);

            % Ensure roiMap is 2D for indexing with central masks
            % If 3D label map, we collapse using max over slices (labels preserved).
            if ndims(roiMap) == 3
                roiMap2D = max(roiMap, [], 3);
            else
                roiMap2D = roiMap;
            end

            % Build 4 traces: average over ROIs whose labels fall under each central mask
            nFrames = size(fluorescenceData, 2);
            fluoroTraces = zeros(nFrames, numel(regions)); % frames x 4

            for k = 1:numel(regions)
                regionKey = regions{k};
                if ~isfield(centralMasks, regionKey)
                    warning("centralMasks missing field %s in %s — skipping region.", regionKey, replicateID);
                    continue;
                end

                mask = logical(centralMasks.(regionKey));
                if ~any(mask(:))
                    warning("Central mask empty (%s) in %s", regionKey, replicateID);
                    continue;
                end

                roiLabels = roiMap2D(mask);
                roiLabels = unique(roiLabels(roiLabels > 0));

                if isempty(roiLabels)
                    warning("No ROI labels found under %s mask in %s", regionKey, replicateID);
                    continue;
                end

                roiIdx = ismember(uniqueROIs, roiLabels);

                if any(roiIdx)
                    % mean across matched ROIs => 1 x nFrames, transpose to nFrames x 1
                    fluoroTraces(:, k) = mean(fluorescenceData(roiIdx, :), 1, "omitnan")';
                else
                    warning("No matching ROIs in fluorescenceData for %s in %s", regionKey, replicateID);
                end
            end

            % Correlation across the 4 region-traces
            corrMatrix = corr(fluoroTraces, "Type", corrType, "Rows", "pairwise");

            % Store replicate-level pairwise entries
            for a = 1:numel(regions)
                for b = (a+1):numel(regions)
                    pairwiseCorrelationData(end+1, :) = { ...
                        replicateID, groupName, labelNames{a}, labelNames{b}, corrMatrix(a,b)};
                end
            end

            % Stack for group-level mean
            if strcmpi(groupName, "Control")
                controlMatrixStack = cat(3, controlMatrixStack, corrMatrix);
            else
                testMatrixStack = cat(3, testMatrixStack, corrMatrix);
            end

        catch ME
            warning("Skipping replicate %s: %s", replicateID, ME.message);
        end
    end
end

%% ============================
% 3) SAVE CSV OUTPUTS
% ============================

pairwiseTable = cell2table(pairwiseCorrelationData, ...
    "VariableNames", {"Replicate","Group","Region1","Region2","Correlation"});
writetable(pairwiseTable, fullfile(outputFolder, "Pairwise_CentralROI_Correlations.csv"));

meanControl = mean(controlMatrixStack, 3, "omitnan");
meanTest    = mean(testMatrixStack, 3, "omitnan");

controlTable = array2table(meanControl, "VariableNames", labelNames, "RowNames", labelNames);
testTable    = array2table(meanTest,    "VariableNames", labelNames, "RowNames", labelNames);

writetable(controlTable, fullfile(outputFolder, "Control_Group_CentralROI_CorrelationMatrix.csv"), "WriteRowNames", true);
writetable(testTable,    fullfile(outputFolder, "Test_Group_CentralROI_CorrelationMatrix.csv"),    "WriteRowNames", true);

% Save separate CSVs for each region pair
uniquePairs = unique(pairwiseTable(:, {"Region1","Region2"}), "rows");
for i = 1:height(uniquePairs)
    region1 = uniquePairs.Region1{i};
    region2 = uniquePairs.Region2{i};
    pairName = sprintf("%s_vs_%s", strrep(region1," ","_"), strrep(region2," ","_"));
    subset = pairwiseTable(strcmp(pairwiseTable.Region1, region1) & strcmp(pairwiseTable.Region2, region2), :);
    writetable(subset, fullfile(outputFolder, pairName + "_Correlation.csv"));
end

%% ============================
% 4) HEATMAP PLOTS (GLOBAL SCALE)
% ============================

rawMin = min([meanControl(:); meanTest(:)]);
globalMin = floor(rawMin * 10) / 10;
globalMax = 1;

plot_corr_heatmap(meanControl, labelNames, ...
    "Control Group — Central ROI Correlation", ...
    fullfile(outputFolder, "Control_Group_CentralROI_CorrelationMatrix.png"), ...
    globalMin, globalMax);

plot_corr_heatmap(meanTest, labelNames, ...
    "Test Group — Central ROI Correlation", ...
    fullfile(outputFolder, "Test_Group_CentralROI_CorrelationMatrix.png"), ...
    globalMin, globalMax);

fprintf("Saved all group-level central ROI correlation outputs to: %s\n", outputFolder);

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

function plot_corr_heatmap(corrMatrix, labels, titleStr, filename, globalMin, globalMax)
    % Colorblind-safe reddish to yellow palette
    cmap = [ ...
        165,0,38;
        215,48,39;
        244,109,67;
        253,174,97;
        254,224,139;
        255,255,191] / 255;

    figure("Visible","off");
    imagesc(corrMatrix);
    colormap(cmap);
    colorbar;
    caxis([globalMin, globalMax]);

    xticks(1:numel(labels)); yticks(1:numel(labels));
    xticklabels(labels); yticklabels(labels);
    xlabel("Region"); ylabel("Region");
    title(titleStr, "FontWeight","bold");

    cmapSize = size(cmap, 1);
    for row = 1:numel(labels)
        for col = 1:numel(labels)
            val = corrMatrix(row, col);

            % protect against NaN
            if isnan(val)
                txt = "NaN";
                textColor = "black";
            else
                colorIndex = round((val - globalMin) / (globalMax - globalMin) * (cmapSize - 1)) + 1;
                colorIndex = max(min(colorIndex, cmapSize), 1);

                rgb = cmap(colorIndex, :);
                brightness = 0.2126 * rgb(1) + 0.7152 * rgb(2) + 0.0722 * rgb(3);
                textColor = "black";
                if brightness < 0.5, textColor = "white"; end
                txt = sprintf("%.2f", val);
            end

            text(col, row, txt, ...
                "HorizontalAlignment","center", ...
                "VerticalAlignment","middle", ...
                "FontSize",10, "FontWeight","bold", ...
                "Color", textColor);
        end
    end

    axis square;
    set(gca, "FontSize", 10);
    saveas(gcf, filename);
    close(gcf);
end
