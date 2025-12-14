%% Step 4 — Whole-brain FFT analysis (power spectrum + band power)
%
% Public/GitHub-ready script (code-only release).
% Uses Step 2 outputs:
%   results/step2_raw_trace_collection/Control_Group_Data.mat  (variable: controlData)
%   results/step2_raw_trace_collection/Test_Group_Data.mat     (variable: testData)
%
% Expected structure in controlData/testData (per replicate):
%   - fluorescenceData : [nROIs x nFrames]
%   - replicateID      : string/char (optional; will be generated if missing)
%
% Outputs:
%   results/step4_fft/
%     - PowerSpectrum_WholeBrain_<Group>_<ReplicateID>.png (per replicate)
%     - Overlay_All_PowerSpectra_Control_vs_Test.png
%     - Overlay_Mean_SD_PowerSpectra_Control_vs_Test.png
%     - Overlay_Mean_95CI_PowerSpectra_Control_vs_Test.png
%     - WholeBrain_FrequencyBand_Summary_with_LogLogPower.csv
%
% Notes:
%   - fps defaults to 2 Hz (edit to match acquisition)
%   - spectra are interpolated onto a fixed frequency grid (500 points)
%   - "band power" is computed via trapezoidal integration over the band
%   - y-axis is log scale; values are protected from log(0) by a minimum

%% ================================
% 0) CLEAN START (optional)
% ================================
clearvars;
close all;
clc;

%% ================================
% 1) USER CONFIGURATION
% ================================

projectRoot = pwd;

% Step 2 input directory (GitHub-safe default)
step2Dir = fullfile(projectRoot, "results", "step2_raw_trace_collection");
controlDataPath = fullfile(step2Dir, "Control_Group_Data.mat");
testDataPath    = fullfile(step2Dir, "Test_Group_Data.mat");

% Output directory for Step 4
outputDir = fullfile(projectRoot, "results", "step4_fft");
if ~exist(outputDir, "dir")
    mkdir(outputDir);
end

% Acquisition
fps = 2;

% Frequency grid (used for interpolation and group overlays)
nFreqPoints    = 500;
allFrequencies = linspace(0, fps/2, nFreqPoints);

% Plot limits
xLimits = [0, 1];
yLimits = [1e-7, 1e-1];

% Color palette (kept as in your original)
colorControl      = [0.0000, 0.4470, 0.7410];  % Blue
colorTest         = [0.8350, 0.3686, 0.0000];  % Orange
colorControlLight = [0.6, 0.8, 1.0];
colorTestLight    = [1.0, 0.7, 0.5];

% Frequency bands (edit/add as needed)
bands = struct();
bands.FullBand = [0.01, 1.0];

% Minimum y to avoid log(0)
minY = 1e-8;

%% ================================
% 2) LOAD STEP 2 DATA
% ================================
disp("Loading control and test group data (Step 2 outputs) for FFT analysis...");

assert(exist(controlDataPath, "file") == 2, "Missing file: %s", controlDataPath);
assert(exist(testDataPath, "file") == 2, "Missing file: %s", testDataPath);

S1 = load(controlDataPath);
S2 = load(testDataPath);

assert(isfield(S1, "controlData"), "Control MAT must contain variable 'controlData'.");
assert(isfield(S2, "testData"),    "Test MAT must contain variable 'testData'.");

controlData = S1.controlData;
testData    = S2.testData;

controlData = ensureReplicateIDs(controlData, "Control");
testData    = ensureReplicateIDs(testData, "Test");

%% ================================
% 3) PER-REPLICATE FFT + BAND POWER
% ================================

nControl = numel(controlData);
nTest    = numel(testData);

summaryData = cell(nControl + nTest, 4); % Group | ReplicateID | FullBandPower | Variance

allPowerSpectra = zeros(nControl + nTest, nFreqPoints);

% ---- Controls ----
for i = 1:nControl
    repID = string(controlData(i).replicateID);
    fluorescence = controlData(i).fluorescenceData;
    [spec, totalBandPower] = computeReplicateSpectrumAndBandPower( ...
        fluorescence, fps, allFrequencies, bands.FullBand);

    % Fix Nyquist artifact (kept from your original)
    spec(end) = spec(end-1);

    allPowerSpectra(i, :) = spec;

    % Save per-replicate spectrum plot
    outFig = fullfile(outputDir, "PowerSpectrum_WholeBrain_Control_" + repID + ".png");
    savePowerSpectrumFigure(allFrequencies, spec, colorControl, xLimits, yLimits, ...
        "Power Spectrum — Whole Brain Control — " + repID, outFig);

    % Replicate variance of mean trace (kept from your original)
    replicateVariance = var(mean(fluorescence, 1, "omitnan"));

    summaryData{i, 1} = "Control";
    summaryData{i, 2} = char(repID);
    summaryData{i, 3} = totalBandPower;
    summaryData{i, 4} = replicateVariance;
end

% ---- Tests ----
for i = 1:nTest
    repID = string(testData(i).replicateID);
    fluorescence = testData(i).fluorescenceData;
    [spec, totalBandPower] = computeReplicateSpectrumAndBandPower( ...
        fluorescence, fps, allFrequencies, bands.FullBand);

    spec(end) = spec(end-1);

    rowIdx = nControl + i;
    allPowerSpectra(rowIdx, :) = spec;

    outFig = fullfile(outputDir, "PowerSpectrum_WholeBrain_Test_" + repID + ".png");
    savePowerSpectrumFigure(allFrequencies, spec, colorTest, xLimits, yLimits, ...
        "Power Spectrum — Whole Brain Test — " + repID, outFig);

    replicateVariance = var(mean(fluorescence, 1, "omitnan"));

    summaryData{rowIdx, 1} = "Test";
    summaryData{rowIdx, 2} = char(repID);
    summaryData{rowIdx, 3} = totalBandPower;
    summaryData{rowIdx, 4} = replicateVariance;
end

%% ================================
% 4) GROUP OVERLAYS (Mean ± SEM, SD, CI)
% ================================

controlSpectra = allPowerSpectra(1:nControl, :);
testSpectra    = allPowerSpectra(nControl+1:end, :);

controlMean = smoothdata(mean(controlSpectra, 1), "movmean", 5);
controlSEM  = std(controlSpectra, 0, 1) / sqrt(max(nControl,1));

testMean = smoothdata(mean(testSpectra, 1), "movmean", 5);
testSEM  = std(testSpectra, 0, 1) / sqrt(max(nTest,1));

% ---- Overlay: all replicates + Mean ± SEM ----
outOverlaySEM = fullfile(outputDir, "Overlay_All_PowerSpectra_Control_vs_Test.png");
plotOverlayAllAndMeanBand(allFrequencies, controlSpectra, testSpectra, ...
    controlMean, controlSEM, testMean, testSEM, ...
    colorControl, colorControlLight, colorTest, colorTestLight, ...
    xLimits, yLimits, outOverlaySEM, "Overlay of Power Spectra — Control vs Test (Mean ± SEM)");

% ---- Mean ± SD ----
controlSD = std(controlSpectra, 0, 1);
testSD    = std(testSpectra, 0, 1);

outOverlaySD = fullfile(outputDir, "Overlay_Mean_SD_PowerSpectra_Control_vs_Test.png");
plotOverlayMeanBand(allFrequencies, controlMean, controlSD, testMean, testSD, ...
    colorControl, colorTest, xLimits, yLimits, minY, outOverlaySD, "Mean ± SD of Power Spectra — Control vs Test");

% ---- Mean ± 95% CI ----
controlCI = controlSEM * 1.96;
testCI    = testSEM * 1.96;

outOverlayCI = fullfile(outputDir, "Overlay_Mean_95CI_PowerSpectra_Control_vs_Test.png");
plotOverlayMeanBand(allFrequencies, controlMean, controlCI, testMean, testCI, ...
    colorControl, colorTest, xLimits, yLimits, minY, outOverlayCI, "Mean ± 95% CI of Power Spectra — Control vs Test");

%% ================================
% 5) SUMMARY CSV
% ================================

tableHeaders = { ...
    "Group", "Replicate", ...
    "FullBand_0_01to1Hz_mV2", ...
    "SignalVariance_DeltaFoverF_squared"};

summaryTable = cell2table(summaryData, "VariableNames", tableHeaders);
writetable(summaryTable, fullfile(outputDir, "WholeBrain_FrequencyBand_Summary_with_LogLogPower.csv"));

disp("Step 4 complete. FFT outputs saved.");

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

function [totalPowerSpectrum, totalBandPower] = computeReplicateSpectrumAndBandPower(fluorescence, fps, freqGrid, band)
    % fluorescence: [nROIs x nFrames]
    nROIs = size(fluorescence, 1);

    totalPowerSpectrum = zeros(1, numel(freqGrid));
    totalBandPower = 0;

    validCount = 0;

    for roiIdx = 1:nROIs
        signal = fluorescence(roiIdx, :);

        % Skip invalid or flat signals (kept from your original)
        if any(isnan(signal)) || any(isinf(signal)) || range(signal) < 1e-6
            continue;
        end

        [frequencies, powerSpectrum] = computeFFT(signal, fps);

        % Sum interpolated spectra onto common grid
        totalPowerSpectrum = totalPowerSpectrum + interp1( ...
            frequencies, powerSpectrum, freqGrid, "linear", "extrap");

        % Sum band power
        totalBandPower = totalBandPower + computeTotalBandPower( ...
            frequencies, powerSpectrum, band, true);

        validCount = validCount + 1;
    end

    % Average spectrum across ROIs (use validCount to avoid bias if many ROIs skipped)
    if validCount > 0
        totalPowerSpectrum = totalPowerSpectrum / validCount;
    else
        totalPowerSpectrum(:) = NaN;
        totalBandPower = NaN;
    end
end

function [frequencies, powerSpectrum_mV2] = computeFFT(signal, fps)
    L = length(signal);
    Y = fft(signal);
    P2 = abs(Y / L);
    P1 = P2(1:floor(L/2)+1);
    if numel(P1) > 2
        P1(2:end-1) = 2 * P1(2:end-1);
    end
    frequencies = fps * (0:floor(L/2)) / L;
    powerSpectrum_mV2 = P1.^2;
end

function bandPower = computeTotalBandPower(frequencies, powerSpectrum, band, isLastBand)
    if isLastBand
        bandIndices = (frequencies >= band(1)) & (frequencies <= band(2));
    else
        bandIndices = (frequencies >= band(1)) & (frequencies < band(2));
    end

    if any(bandIndices) && numel(frequencies(bandIndices)) > 1
        bandPower = trapz(frequencies(bandIndices), powerSpectrum(bandIndices));
    else
        bandPower = 0;
    end
end

function savePowerSpectrumFigure(freq, spec, lineColor, xLimits, yLimits, figTitle, outFile)
    fig = figure("Visible", "off");
    semilogy(freq, spec, "Color", lineColor, "LineWidth", 2);

    title(figTitle);
    xlabel("Frequency (Hz)");
    ylabel("Power (\DeltaF/F)^2 (Log Scale)");

    xlim(xLimits);
    ylim(yLimits);

    ax = gca;
    set(ax, ...
        "TickDir", "out", ...
        "TickLength", [0.015 0.015], ...
        "YMinorTick", "off", ...
        "XMinorTick", "off", ...
        "YScale", "log", ...
        "XColor", "k", "YColor", "k", ...
        "Color", "w", ...
        "XGrid", "off", "YGrid", "off", ...
        "Box", "off");

    set(gcf, "Color", "w");

    saveas(fig, outFile);
    close(fig);
end

function plotOverlayAllAndMeanBand(freq, controlSpectra, testSpectra, controlMean, controlSEM, testMean, testSEM, ...
    colorControl, colorControlLight, colorTest, colorTestLight, xLimits, yLimits, outFile, figTitle)

    fig = figure("Visible", "off");
    hold on;

    % all replicate lines
    for i = 1:size(controlSpectra, 1)
        semilogy(freq, controlSpectra(i, :), "Color", colorControlLight, "LineWidth", 1, "HandleVisibility", "off");
    end
    for i = 1:size(testSpectra, 1)
        semilogy(freq, testSpectra(i, :), "Color", colorTestLight, "LineWidth", 1, "HandleVisibility", "off");
    end

    % Mean ± SEM shading (protect against <=0 on log scale)
    ctrl_upper = max(controlMean + controlSEM, 1e-8);
    ctrl_lower = max(controlMean - controlSEM, 1e-8);
    test_upper = max(testMean + testSEM, 1e-8);
    test_lower = max(testMean - testSEM, 1e-8);

    fill([freq, fliplr(freq)], [ctrl_upper, fliplr(ctrl_lower)], colorControl, ...
        "FaceAlpha", 0.2, "EdgeColor", "none", "HandleVisibility", "off");
    h1 = semilogy(freq, controlMean, "Color", colorControl, "LineWidth", 3, "DisplayName", "Control Mean ± SEM");

    fill([freq, fliplr(freq)], [test_upper, fliplr(test_lower)], colorTest, ...
        "FaceAlpha", 0.2, "EdgeColor", "none", "HandleVisibility", "off");
    h2 = semilogy(freq, testMean, "Color", colorTest, "LineWidth", 3, "DisplayName", "Test Mean ± SEM");

    xlabel("Frequency (Hz)");
    ylabel("Power (\DeltaF/F)^2 (Log Scale)");
    title(figTitle);

    ax = gca;
    set(ax, ...
        "TickDir", "out", ...
        "TickLength", [0.015 0.015], ...
        "YMinorTick", "off", ...
        "XMinorTick", "off", ...
        "YScale", "log", ...
        "XColor", "k", "YColor", "k", ...
        "Color", "w", ...
        "XGrid", "off", "YGrid", "off", ...
        "Box", "off");
    set(gcf, "Color", "w");

    xlim(xLimits);
    ylim(yLimits);

    legend([h1, h2], "Location", "northeast");
    hold off;

    saveas(fig, outFile);
    close(fig);
end

function plotOverlayMeanBand(freq, mean1, bandWidth1, mean2, bandWidth2, color1, color2, xLimits, yLimits, minY, outFile, figTitle)
    fig = figure("Visible", "off");
    hold on;

    upper1 = max(mean1 + bandWidth1, minY);
    lower1 = max(mean1 - bandWidth1, minY);
    upper2 = max(mean2 + bandWidth2, minY);
    lower2 = max(mean2 - bandWidth2, minY);

    fill([freq, fliplr(freq)], [upper1, fliplr(lower1)], color1, "FaceAlpha", 0.3, "EdgeColor", "none");
    fill([freq, fliplr(freq)], [upper2, fliplr(lower2)], color2, "FaceAlpha", 0.3, "EdgeColor", "none");

    h1 = plot(freq, mean1, "-", "Color", color1, "LineWidth", 2.5);
    h2 = plot(freq, mean2, "-", "Color", color2, "LineWidth", 2.5);

    set(gca, "YScale", "log");
    xlabel("Frequency (Hz)");
    ylabel("Power (\DeltaF/F)^2 (Log Scale)");
    title(figTitle);

    ax = gca;
    ax.TickDir = "out";
    ax.TickLength = [0.015 0.015];
    ax.YMinorTick = "off";
    ax.XMinorTick = "off";
    ax.XColor = "k";
    ax.YColor = "k";
    ax.Color = "w";
    ax.XGrid = "off";
    ax.YGrid = "off";
    ax.Box = "off";

    xlim(xLimits);
    ylim(yLimits);
    legend([h1, h2], {"Control", "Test"}, "Location", "northeast");
    set(gcf, "Color", "w");

    saveas(fig, outFile);
    close(fig);
end
