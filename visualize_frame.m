function visualize_frame(txData, chData, rxData, cfg)
% =========================================================================
% VISUALIZE_FRAME — Plot TX/RX waveform and constellations for one frame
% =========================================================================
% Produces three figures:
%
%   Figure 1 — Waveform (2x2):
%       TX I channel | RX I channel
%       TX Q channel | RX Q channel
%
%   Figure 2 — Constellation (1x3):
%       TX | RX before LS equalization | RX after LS equalization
%
%   Figure 3 — LS channel estimates per pilot block (1x2):
%       Re(h_LS) vs pilot block index | Im(h_LS) vs pilot block index
%
% Controlled by cfg fields:
%   cfg.Visualize      = true/false   — master switch
%   cfg.VisualizeFrame = 'random'     — pick random frame
%                      = N (integer)  — pick specific frame index
%
% All figures are saved as PNG to cfg.OutputDir.
% =========================================================================

fprintf('    [VIZ] Generating plots...\n');

%% -----------------------------------------------------------------------
%  Select frame index
% -----------------------------------------------------------------------
if ischar(cfg.VisualizeFrame) && strcmpi(cfg.VisualizeFrame,'random')
    frameIdx = randi(cfg.NumFrames);
else
    frameIdx = max(1, min(round(cfg.VisualizeFrame), cfg.NumFrames));
end
fprintf('    [VIZ] Frame selected: %d / %d\n', frameIdx, cfg.NumFrames);

%% -----------------------------------------------------------------------
%  Extract raw waveform samples for this frame
% -----------------------------------------------------------------------
totalSamples    = length(txData.waveform);
samplesPerFrame = floor(totalSamples / cfg.NumFrames);
fStart          = (frameIdx-1)*samplesPerFrame + 1;
fEnd            = min(fStart + samplesPerFrame - 1, totalSamples);

txFrame   = txData.waveform(fStart:fEnd);
rxFrame   = chData.rxWaveform(fStart:fEnd);
sampleIdx = fStart:fEnd;

%% -----------------------------------------------------------------------
%  Extract symbols (downsample) for constellation
% -----------------------------------------------------------------------
mkFilter = @() comm.RaisedCosineReceiveFilter(...
    'RolloffFactor',         cfg.RolloffFactor, ...
    'FilterSpanInSymbols',   cfg.FilterSpanInSymbols, ...
    'InputSamplesPerSymbol', cfg.SamplesPerSymbol, ...
    'DecimationFactor',      cfg.SamplesPerSymbol);

delay     = cfg.FilterSpanInSymbols / 2;

f1 = mkFilter();
txSymAll = f1(txData.waveform);
txSymAll = txSymAll(delay+1:end);
release(f1);

f2 = mkFilter();
rxSymAll = f2(chData.rxWaveform);
rxSymAll = rxSymAll(delay+1:end);
release(f2);

minLen   = min(length(txSymAll), length(rxSymAll));
txSymAll = txSymAll(1:minLen);
rxSymAll = rxSymAll(1:minLen);

symPerFrame = floor(minLen / cfg.NumFrames);
sStart      = (frameIdx-1)*symPerFrame + 1;
sEnd        = min(sStart + symPerFrame - 1, minLen);

txSym   = txSymAll(sStart:sEnd);
rxSym   = rxSymAll(sStart:sEnd);
rxSymEq = rxSym / rxData.h_LS_frame;    % LS equalized

%% -----------------------------------------------------------------------
%  FIGURE 1: TX and RX Waveforms
% -----------------------------------------------------------------------
fig1 = figure('Name', sprintf('Waveform | %s | Frame %d', ...
    cfg.CurrentMODCODName, frameIdx), ...
    'Position', [50 300 1200 480], 'Color', 'w');

% TX I
subplot(2,2,1);
plot(sampleIdx, real(txFrame), 'b', 'LineWidth', 0.9);
grid on; box on;
xlabel('Sample Index'); ylabel('Amplitude');
title(sprintf('TX — I channel  (Frame %d)', frameIdx));
xlim([sampleIdx(1) sampleIdx(end)]);

% TX Q
subplot(2,2,3);
plot(sampleIdx, imag(txFrame), 'b', 'LineWidth', 0.9);
grid on; box on;
xlabel('Sample Index'); ylabel('Amplitude');
title('TX — Q channel');
xlim([sampleIdx(1) sampleIdx(end)]);

% RX I
subplot(2,2,2);
plot(sampleIdx, real(rxFrame), 'r', 'LineWidth', 0.9);
grid on; box on;
xlabel('Sample Index'); ylabel('Amplitude');
title(sprintf('RX — I channel  (SNR=%.0f dB  |  A_p=%.2f dB)', ...
    chData.snr_dB, chData.rainAtten_dB));
xlim([sampleIdx(1) sampleIdx(end)]);

% RX Q
subplot(2,2,4);
plot(sampleIdx, imag(rxFrame), 'r', 'LineWidth', 0.9);
grid on; box on;
xlabel('Sample Index'); ylabel('Amplitude');
title('RX — Q channel');
xlim([sampleIdx(1) sampleIdx(end)]);

sgtitle(sprintf('DVB-S2 Waveform  |  MODCOD: %s  |  |h_{tilde}| = %.4f  |  Frame: %d/%d', ...
    cfg.CurrentMODCODName, abs(chData.h_tilde), frameIdx, cfg.NumFrames), ...
    'FontSize', 12, 'FontWeight', 'bold');

saveas(fig1, fullfile(cfg.OutputDir, ...
    sprintf('fig1_waveform_%s_fr%d.png', cfg.CurrentMODCODName, frameIdx)));

%% -----------------------------------------------------------------------
%  FIGURE 2: Constellation Diagrams
% -----------------------------------------------------------------------
fig2 = figure('Name', sprintf('Constellation | %s | Frame %d', ...
    cfg.CurrentMODCODName, frameIdx), ...
    'Position', [80 80 1300 420], 'Color', 'w');

% Reference QPSK points
ref = (1/sqrt(2)) * [1+1j, 1-1j, -1+1j, -1-1j];

% --- TX ---
subplot(1,3,1);
scatter(real(txSym), imag(txSym), 10, [0.2 0.4 0.8], 'filled', ...
    'MarkerFaceAlpha', 0.35);
hold on;
scatter(real(ref), imag(ref), 150, 'k', 'x', 'LineWidth', 2.5);
hold off;
grid on; axis equal; box on;
xlabel('I'); ylabel('Q');
title('TX Constellation');
xlim([-2 2]); ylim([-2 2]);
legend('TX','Reference','Location','se','FontSize',9);

% --- RX before equalization ---
subplot(1,3,2);
rxLim = max(1.5, max(abs([real(rxSym); imag(rxSym)])) * 1.15);
scatter(real(rxSym), imag(rxSym), 10, [0.8 0.2 0.2], 'filled', ...
    'MarkerFaceAlpha', 0.35);
hold on;
scatter(real(ref), imag(ref), 150, 'k', 'x', 'LineWidth', 2.5);
hold off;
grid on; axis equal; box on;
xlabel('I'); ylabel('Q');
title({'RX Constellation', 'Before LS Equalization'});
xlim([-rxLim rxLim]); ylim([-rxLim rxLim]);
legend('RX','Reference','Location','se','FontSize',9);

% --- RX after LS equalization ---
subplot(1,3,3);
scatter(real(rxSymEq), imag(rxSymEq), 10, [0.1 0.6 0.1], 'filled', ...
    'MarkerFaceAlpha', 0.35);
hold on;
scatter(real(ref), imag(ref), 150, 'k', 'x', 'LineWidth', 2.5);
hold off;
grid on; axis equal; box on;
xlabel('I'); ylabel('Q');
title({'RX Constellation', 'After LS Equalization'});
xlim([-2 2]); ylim([-2 2]);
legend('Equalized','Reference','Location','se','FontSize',9);

sgtitle(sprintf('Constellation  |  MODCOD: %s  |  SNR: %.0f dB  |  A_p: %.2f dB  |  BER(LS): %.5f', ...
    cfg.CurrentMODCODName, chData.snr_dB, chData.rainAtten_dB, rxData.BER_LS), ...
    'FontSize', 12, 'FontWeight', 'bold');

saveas(fig2, fullfile(cfg.OutputDir, ...
    sprintf('fig2_constellation_%s_fr%d.png', cfg.CurrentMODCODName, frameIdx)));

%% -----------------------------------------------------------------------
%  FIGURE 3: LS Estimates per Pilot Block
% -----------------------------------------------------------------------
fig3 = figure('Name', sprintf('LS Estimates | %s | Frame %d', ...
    cfg.CurrentMODCODName, frameIdx), ...
    'Position', [120 150 900 380], 'Color', 'w');

pb = 1:rxData.NumPilotBlocks;

subplot(1,2,1);
plot(pb, real(rxData.h_LS_perBlock), 'b-o', ...
    'LineWidth', 1.4, 'MarkerSize', 4, 'DisplayName', 'Re(h_{LS})');
hold on;
yline(real(chData.h_tilde), 'r--', 'LineWidth', 1.8, ...
    'DisplayName', 'True Re(h_{tilde})');
hold off;
grid on; box on;
xlabel('Pilot Block Index'); ylabel('Value');
title('Re(h_{LS}) per Pilot Block');
legend('Location','best','FontSize',9);

subplot(1,2,2);
plot(pb, imag(rxData.h_LS_perBlock), 'r-o', ...
    'LineWidth', 1.4, 'MarkerSize', 4, 'DisplayName', 'Im(h_{LS})');
hold on;
yline(imag(chData.h_tilde), 'b--', 'LineWidth', 1.8, ...
    'DisplayName', 'True Im(h_{tilde})');
hold off;
grid on; box on;
xlabel('Pilot Block Index'); ylabel('Value');
title('Im(h_{LS}) per Pilot Block');
legend('Location','best','FontSize',9);

sgtitle(sprintf('LS Channel Estimates  |  MODCOD: %s  |  Frame: %d  |  NMSE: %.4e', ...
    cfg.CurrentMODCODName, frameIdx, rxData.NMSE_LS), ...
    'FontSize', 12, 'FontWeight', 'bold');

saveas(fig3, fullfile(cfg.OutputDir, ...
    sprintf('fig3_ls_estimates_%s_fr%d.png', cfg.CurrentMODCODName, frameIdx)));

%% Summary line
fprintf('    [VIZ] Frame %d | |h_tilde|=%.4f | A_p=%.2f dB | BER=%.5f | NMSE=%.4e\n', ...
    frameIdx, abs(chData.h_tilde), chData.rainAtten_dB, rxData.BER_LS, rxData.NMSE_LS);
fprintf('    [VIZ] 3 figures saved to: %s\n', cfg.OutputDir);
end
