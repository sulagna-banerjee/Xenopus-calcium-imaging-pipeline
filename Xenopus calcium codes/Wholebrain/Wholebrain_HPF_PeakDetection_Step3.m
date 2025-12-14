%% Step 3 — HPF-based peak detection + overlay plots + summary tables
%
% Public/GitHub-ready script (code-only release).
% Expects Step 2 outputs containing:
%   - controlData struct array with fields: fluorescenceData, replicateID (recommended)
%   - testData struct array with fields: fluorescenceData, replicateID (recommended)
%
% Folder layout (suggested):
% project_root/
% ├── scripts/
% │   └── step3_hpf_peak_detection_and_plots.m
% ├── results/
% │   ├── step2_raw_trace_collection/
% │   │   ├── Control_Group_Data.mat
% │   │   └── Test_Group_Data.mat
% │   └── step3_hpf_peaks/
%
% Outputs:
%   results/step3_hpf_peaks/
%     - Significant_Peaks_Summary_with_HPF.csv
%     - Individual_Peak_Details.csv
%     - <Group>_<ReplicateID>_FinalOverlayPlot.png (per replicate)
%     - Standalone_Legend_AllSignals.png/.pdf
%
% Notes:
%   - fps defaults to 2 (edit to match your acquisition rate)
%   - High-pass filter cutoff fixed at 0.005 Hz

%% ============================
% 0) USER CONFIGURATION
% ============================

projectRoot = pwd;

% ---- Input Step 2 outputs ----
step2Dir = fullfile(projectRoot, "results", "step2_raw_trace_collection");
controlMatPath = fullfile(step2Dir, "Control_Group_Data.mat");
testMatPath    = fullfile(step2Dir, "Test_Group_Data.mat");

% ---- Output directory for Step 3 ----
outputFolder = fullfile(projectRoot, "results", "step3_hpf_peaks");
if ~exist(outputFolder, "dir")
    mkdir(outputFolder);
end

% ---- Acquisition rate ----
fps = 2;

% ---- Cutoff configuration ----
% Choose which group sets the cutoff distribution: "control" or "test"
cutoffChoice = "control";     % <-- change to "test" if needed
sdMultiplier = 3;             % e.g., 3 × SD

% ---- Plot scaling ----
globalYMin = -10; % fixed
% globalYMax computed automatically from data

%% ============================
% 1) LOAD DATA
% ============================

assert(exist(controlMatPath, "file") == 2, "Missing file: %s", controlMatPath);
assert(exist(testMatPath, "file") == 2, "Missing file: %s", testMatPath);

controlStruct = load(controlMatPath);
testStruct    = load(testMatPath);

if isfield(controlStruct, "controlData")
    controlData = controlStruct.controlData;
else
    error("Control MAT does not contain 'controlData'.");
end

if isfield(testStruct, "testData")
    testData = testStruct.testData;
else
    error("Test MAT does not contain 'testData'.");
end

disp("Loaded Step 2 data successfully.");

% Ensure replicateID exists; if not, create safe placeholders
controlData = ensureReplicateIDs(controlData, "Control");
testData    = ensureReplicateIDs(testData, "Test");

%% ============================
% 2) BUILD DISTRIBUTION FOR GLOBAL CUTOFF
% ============================

if strcmpi(cutoffChoice, "control")
    allFilteredSignals = collectFilteredMeans(controlData, fps);
    baseGroup = "Control";
elseif strcmpi(cutoffChoice, "test")
    allFilteredSignals = collectFilteredMeans(testData, fps);
    baseGroup = "Test";
else
    error('cutoffChoice must be "control" or "test".');
end

globalCutoff = sdMultiplier * std(allFilteredSignals, "omitnan");
disp("Global cutoff computed from " + baseGroup + " group using " + sdMultiplier + " × SD.");
disp("Cutoff value: " + num2str(globalCutoff, "%.3f"));

%% ============================
% 3) COMPUTE GLOBAL Y-MAX (CONSISTENT PLOTS)
% ============================

allData = [controlData, testData];
globalYMax = computeGlobalYMax(allData, fps) + 1;

%% ============================
% 4) SUMMARISE PEAKS + MAKE OVERLAY PLOTS
% ============================

[controlSummary, controlPeaks] = summarizeSignificantPeaks( ...
    controlData, "Control", fps, globalCutoff, outputFolder, globalYMin, globalYMax);

[testSummary, testPeaks] = summarizeSignificantPeaks( ...
    testData, "Test", fps, globalCutoff, outputFolder, globalYMin, globalYMax);

combinedSummary = [controlSummary; testSummary];
writetable(combinedSummary, fullfile(outputFolder, "Significant_Peaks_Summary_with_HPF.csv"));

peakDetails = [controlPeaks; testPeaks];
writetable(peakDetails, fullfile(outputFolder, "Individual_Peak_Details.csv"));

disp("Summary and individual peak details saved.");

%% ============================
% 5) STANDALONE LEGEND EXPORT
% ============================

exportStandaloneLegend(outputFolder);
disp("Standalone legend saved.");

disp("Step 3 complete.");

%% ============================================================
% FUNCTIONS
%% ============================================================

function data = ensureReplicateIDs(data, groupLabel)
    for i = 1:numel(data)
        if ~isfield(data(i), "replicateID") || strlength(string(data(i).replicateID)) == 0
            data(i).replicateID = groupLabel + "_" + string(i);
        end
    end
end

function allFilteredSignals = collectFilteredMeans(data, fps)
    allFilteredSignals = [];
    for i = 1:numel(data)
        sig = mean(data(i).fluorescenceData, 1, "omitnan");
        filt = highPassFilter(sig, fps);
        allFilteredSignals = [allFilteredSignals, filt];
    end
end

function globalYMax = computeGlobalYMax(allData, fps)
    globalYMax = -inf;
    for i = 1:numel(allData)
        sig  = mean(allData(i).fluorescenceData, 1, "omitnan");
        filt = highPassFilter(sig, fps);
        globalYMax = max(globalYMax, max([sig, filt], [], "omitnan"));
    end
end

function [replicateSummary, peakDetails] = summarizeSignificantPeaks(data, groupLabel, fps, cutoff, outputFolder, yMin, yMax)

    replicateSummary = table;
    peakDetails = table;

    for i = 1:numel(data)
        replicateID = string(data(i).replicateID);

        signal = mean(data(i).fluorescenceData, 1, "omitnan");
        filteredSig = highPassFilter(signal, fps);

        signal = signal(:);
        filteredSig = filteredSig(:);
        time = (1:length(signal)) / fps;

        [pks, locs, w, prom] = findpeaks(filteredSig, ...
            "MinPeakProminence", std(filteredSig) * 0.5, ...
            "MinPeakWidth", round(fps * 5), ...
            "MinPeakDistance", round(fps * 5));

        significantIdx = pks > cutoff;

        % Summary metrics
        sigCount        = sum(significantIdx);
        meanSigAmp      = mean(pks(significantIdx), "omitnan");
        meanSigWidth    = mean(w(significantIdx) / fps, "omitnan");
        meanProminence  = mean(prom(significantIdx), "omitnan");
        hfp             = numel(pks) / (length(signal) / fps);

        replicateSummary = [replicateSummary; table( ...
            {groupLabel}, {char(replicateID)}, ...
            sigCount, meanSigAmp, hfp, meanSigWidth, meanProminence, ...
            "VariableNames", {"Group","ReplicateID","SignificantPeaks", ...
                              "MeanSignificantAmplitude","HFP","MeanPeakWidthSec","MeanProminence"})];

        % Peak details table
        if any(significantIdx)
            times = (locs(significantIdx) / fps);   times = times(:);
            amps  = pks(significantIdx);            amps = amps(:);
            widths = (w(significantIdx) / fps);     widths = widths(:);
            prominences = prom(significantIdx);     prominences = prominences(:);

            n = numel(times);
            if all([numel(amps), numel(widths), numel(prominences)] == n)
                tempDetails = table( ...
                    repmat({groupLabel}, n, 1), ...
                    repmat({char(replicateID)}, n, 1), ...
                    times, amps, widths, prominences, ...
                    "VariableNames", {"Group","ReplicateID","TimeSec","Amplitude","WidthSec","Prominence"});

                peakDetails = [peakDetails; tempDetails];
            else
                warning("Size mismatch in peak arrays for %s — skipping peak details.", replicateID);
            end
        end

        % Overlay Plot (significant peaks marked with black arrows)
        peakTimes = locs(significantIdx) / fps;
        adjustedPeakY = pks(significantIdx) + 1;

        fig = figure("Units", "normalized", "Position", [0.1, 0.1, 0.9, 0.85], "Visible", "off");
        hold on;

        plot(time, signal, ":", "Color", [0.5 0.5 0.5], "LineWidth", 1.0);

        if strcmpi(groupLabel, "Control")
            plot(time, filteredSig, "Color", [1, 0.4, 0], "LineWidth", 1.5);
        else
            plot(time, filteredSig, "Color", [0.4 0 0.6], "LineWidth", 1.5);
        end

        yline(cutoff, "-", "Color", [0 0 0], "LineWidth", 2);
        scatter(peakTimes, adjustedPeakY, 90, "k", "v", "filled", "HandleVisibility", "off");

        xlabel("Time (seconds)", "FontSize", 24, "FontWeight", "bold");
        ylabel("Fluorescence (\DeltaF/F%)", "FontSize", 24, "FontWeight", "bold");
        title("Raw + Filtered Signal — " + groupLabel + " " + replicateID, "FontSize", 26, "FontWeight", "bold");

        ylim([yMin, yMax]);
        xlim([0, time(end)]);
        set(gca, "FontSize", 18, "FontWeight", "bold", "TickDir", "out");
        box off;

        outFile = fullfile(outputFolder, groupLabel + "_" + replicateID + "_FinalOverlayPlot.png");
        saveas(fig, outFile);
        close(fig);
    end
end

function filteredSig = highPassFilter(signal, fps)
    cutoffHz = 0.005;
    [b, a] = butter(2, cutoffHz / (fps / 2), "high");
    filteredSig = filtfilt(b, a, signal);
end

function exportStandaloneLegend(outputFolder)
    fig = figure("Color", "w", "Units", "inches", "Position", [1, 1, 3, 3]);
    axes("Position", [0 0 1 1]); axis off; hold on;

    hRaw      = plot(nan, nan, ":", "Color", [0.5 0.5 0.5], "LineWidth", 1.2);
    hFiltCtrl = plot(nan, nan, "-", "Color", [1 0.4 0], "LineWidth", 1.8);
    hFiltTest = plot(nan, nan, "-", "Color", [0.4 0 0.6], "LineWidth", 1.8);
    hCutoff   = plot(nan, nan, "-", "Color", [0 0 0], "LineWidth", 2);
    hPeak     = plot(nan, nan, "v", "MarkerEdgeColor", "k", ...
        "MarkerFaceColor", "k", "MarkerSize", 10);

    legend([hRaw, hFiltCtrl, hFiltTest, hCutoff, hPeak], ...
        {"Raw Signal", "Control Filtered Signal", "CRISPant Filtered Signal", "Global Cutoff", "Significant Peaks"}, ...
        "Box", "off", "FontSize", 14, "Location", "southoutside", "Orientation", "vertical");

    saveas(fig, fullfile(outputFolder, "Standalone_Legend_AllSignals.png"));
    saveas(fig, fullfile(outputFolder, "Standalone_Legend_AllSignals.pdf"));
    close(fig);
end
