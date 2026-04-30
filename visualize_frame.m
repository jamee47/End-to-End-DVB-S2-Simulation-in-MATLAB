function visualize_frame(txData, chData, rxData, cfg)
% =========================================================================
% VISUALIZE_FRAME — Waveform and constellation plots for one DVB-S2 frame
% =========================================================================
% Produces three figures for a single selected frame:
%
%   Figure 1 — Waveform (2 rows x 2 cols):
%       Top-left  : TX I (in-phase)     Top-right  : RX I (in-phase)
%       Bot-left  : TX Q (quadrature)   Bot-right  : RX Q (quadrature)
%
%   Figure 2 — Constellation (1 row x 3 cols):
%       TX symbols | RX before equalization | RX after LS equalization
%
%   Figure 3 — LS channel estimates across pilot blocks (1 row x 2 cols):
%       Re(h_LS) vs pilot block | Im(h_LS) vs pilot block
%
% Configuration:
%   cfg.Visualize      = true / false
%   cfg.VisualizeFrame = 'random'  OR  N (integer frame index)
%
% All figures saved as PNG to cfg.OutputDir.
% =========================================================================

fprintf('    [VIZ] Generating single-frame plots...\n');

%% -----------------------------------------------------------------------
%  Step 1: Select which frame to visualize
% -----------------------------------------------------------------------
if ischar(cfg.VisualizeFrame) && strcmpi(cfg.VisualizeFrame, 'random')
    frameIdx = randi(cfg.NumFrames);
else
    frameIdx = max(1, min(round(cfg.VisualizeFrame), cfg.NumFrames));
end
fprintf('    [VIZ] Frame: %d / %d\n', frameIdx, cfg.NumFrames);

%% -----------------------------------------------------------------------
%  Step 2: Slice exactly one frame from the raw waveform
%  The waveform is at cfg.SamplesPerSymbol samples/symbol.
%  We slice raw samples — no filtering at this stage.
% -----------------------------------------------------------------------
SPS             = cfg.SamplesPerSymbol;
totalSamples    = length(txData.waveform);
samplesPerFrame = floor(totalSamples / cfg.NumFrames);

fStart = (frameIdx - 1) * samplesPerFrame + 1;
fEnd   = min(fStart + samplesPerFrame - 1, totalSamples);

txFrame   = txData.waveform(fStart:fEnd);   % one frame, SPS samples/symbol
rxFrame   = chData.rxWaveform(fStart:fEnd);
nSamples  = length(txFrame);
tAxis     = (0:nSamples-1);                 % relative sample index

%% -----------------------------------------------------------------------
%  Step 3: Downsample ONE frame to symbol rate for constellation
%  Apply matched filter to the single frame slice only.
%  We feed only the frame samples — not the full waveform.
%  Prepend zeros to handle filter startup (length = filter span * SPS).
% -----------------------------------------------------------------------
filterSpanSamples = cfg.FilterSpanInSymbols * SPS;
pad               = zeros(filterSpanSamples, 1);   % startup padding

mkRxFilter = @() comm.RaisedCosineReceiveFilter( ...
    'RolloffFactor',         cfg.RolloffFactor, ...
    'FilterSpanInSymbols',   cfg.FilterSpanInSymbols, ...
    'InputSamplesPerSymbol', SPS, ...
    'DecimationFactor',      SPS);

% TX frame symbols
f_tx        = mkRxFilter();
txSym_pad   = f_tx([pad; txFrame]);
txSym       = txSym_pad(cfg.FilterSpanInSymbols+1:end);  % remove pad delay
release(f_tx);

% RX frame symbols (before equalization)
f_rx        = mkRxFilter();
rxSym_pad   = f_rx([pad; rxFrame]);
rxSym       = rxSym_pad(cfg.FilterSpanInSymbols+1:end);
release(f_rx);

% RX after LS equalization
rxSymEq = rxSym / rxData.h_LS_frame;

% Trim to same length
minSym  = min([length(txSym), length(rxSym), length(rxSymEq)]);
txSym   = txSym(1:minSym);
rxSym   = rxSym(1:minSym);
rxSymEq = rxSymEq(1:minSym);

%% -----------------------------------------------------------------------
%  FIGURE 1: TX and RX Waveform — I and Q channels
% -----------------------------------------------------------------------
fig1 = figure( ...
    'Name',     sprintf('Waveform | %s | Frame %d', cfg.CurrentMODCODName, frameIdx), ...
    'Position', [50 300 1300 500], ...
    'Color',    'w');

% TX I
subplot(2,2,1);
plot(tAxis, real(txFrame), 'Color', [0.18 0.38 0.72], 'LineWidth', 0.9);
grid on; box on;
xlabel('Sample (relative)');
ylabel('Amplitude');
title(sprintf('TX  —  I channel  (Frame %d)', frameIdx), 'FontSize', 11);
xlim([0 nSamples-1]);

% TX Q
subplot(2,2,3);
plot(tAxis, imag(txFrame), 'Color', [0.18 0.38 0.72], 'LineWidth', 0.9);
grid on; box on;
xlabel('Sample (relative)');
ylabel('Amplitude');
title('TX  —  Q channel', 'FontSize', 11);
xlim([0 nSamples-1]);

% RX I
subplot(2,2,2);
plot(tAxis, real(rxFrame), 'Color', [0.80 0.18 0.18], 'LineWidth', 0.9);
grid on; box on;
xlabel('Sample (relative)');
ylabel('Amplitude');
title(sprintf('RX  —  I channel  |  SNR = %.0f dB  |  A_p = %.2f dB', ...
    chData.snr_dB, chData.rainAtten_dB), 'FontSize', 11);
xlim([0 nSamples-1]);

% RX Q
subplot(2,2,4);
plot(tAxis, imag(rxFrame), 'Color', [0.80 0.18 0.18], 'LineWidth', 0.9);
grid on; box on;
xlabel('Sample (relative)');
ylabel('Amplitude');
title('RX  —  Q channel', 'FontSize', 11);
xlim([0 nSamples-1]);

sgtitle(sprintf( ...
    'DVB-S2 Waveform  |  MODCOD: %s  |  |h_{tilde}| = %.4f  |  Frame %d / %d', ...
    cfg.CurrentMODCODName, abs(chData.h_tilde), frameIdx, cfg.NumFrames), ...
    'FontSize', 12, 'FontWeight', 'bold');

f1name = fullfile(cfg.OutputDir, ...
    sprintf('fig1_waveform_%s_fr%d.png', cfg.CurrentMODCODName, frameIdx));
saveas(fig1, f1name);
fprintf('    [VIZ] Saved: %s\n', f1name);

%% -----------------------------------------------------------------------
%  FIGURE 2: Constellation Diagrams
%  TX  |  RX before EQ  |  RX after LS EQ
% -----------------------------------------------------------------------
fig2 = figure( ...
    'Name',     sprintf('Constellation | %s | Frame %d', cfg.CurrentMODCODName, frameIdx), ...
    'Position', [80 80 1350 430], ...
    'Color',    'w');

% QPSK ideal reference points (DVB-S2 Gray-coded)
refPts = (1/sqrt(2)) * [1+1j, -1+1j, 1-1j, -1-1j];

% --- TX constellation ---
subplot(1,3,1);
scatter(real(txSym), imag(txSym), 12, ...
    repmat([0.18 0.38 0.72], minSym, 1), 'filled', 'MarkerFaceAlpha', 0.4);
hold on;
scatter(real(refPts), imag(refPts), 160, 'k', 'x', 'LineWidth', 2.5);
hold off;
grid on; axis equal; box on;
xlabel('In-phase (I)', 'FontSize', 10);
ylabel('Quadrature (Q)', 'FontSize', 10);
title('TX  —  Clean', 'FontSize', 11);
xlim([-2 2]); ylim([-2 2]);
legend('TX symbols', 'QPSK reference', 'Location', 'se', 'FontSize', 8);

% --- RX before equalization ---
subplot(1,3,2);
rxLim = max(1.8, max(abs([real(rxSym); imag(rxSym)])) * 1.2);
scatter(real(rxSym), imag(rxSym), 12, ...
    repmat([0.80 0.18 0.18], minSym, 1), 'filled', 'MarkerFaceAlpha', 0.4);
hold on;
scatter(real(refPts), imag(refPts), 160, 'k', 'x', 'LineWidth', 2.5);
hold off;
grid on; axis equal; box on;
xlabel('In-phase (I)', 'FontSize', 10);
ylabel('Quadrature (Q)', 'FontSize', 10);
title({'RX  —  Before LS Equalization', ...
    sprintf('(rotated/scaled by h_{tilde})')}, 'FontSize', 11);
xlim([-rxLim rxLim]); ylim([-rxLim rxLim]);
legend('RX symbols', 'QPSK reference', 'Location', 'se', 'FontSize', 8);

% --- RX after LS equalization ---
subplot(1,3,3);
scatter(real(rxSymEq), imag(rxSymEq), 12, ...
    repmat([0.10 0.58 0.18], minSym, 1), 'filled', 'MarkerFaceAlpha', 0.4);
hold on;
scatter(real(refPts), imag(refPts), 160, 'k', 'x', 'LineWidth', 2.5);
hold off;
grid on; axis equal; box on;
xlabel('In-phase (I)', 'FontSize', 10);
ylabel('Quadrature (Q)', 'FontSize', 10);
title({'RX  —  After LS Equalization', ...
    sprintf('(h_{LS} = %.4f + j%.4f)', real(rxData.h_LS_frame), imag(rxData.h_LS_frame))}, ...
    'FontSize', 11);
xlim([-2 2]); ylim([-2 2]);
legend('Equalized', 'QPSK reference', 'Location', 'se', 'FontSize', 8);

sgtitle(sprintf( ...
    'Constellation  |  MODCOD: %s  |  SNR: %.0f dB  |  A_p: %.2f dB  |  BER(LS): %.5f', ...
    cfg.CurrentMODCODName, chData.snr_dB, chData.rainAtten_dB, rxData.BER_LS), ...
    'FontSize', 12, 'FontWeight', 'bold');

f2name = fullfile(cfg.OutputDir, ...
    sprintf('fig2_constellation_%s_fr%d.png', cfg.CurrentMODCODName, frameIdx));
saveas(fig2, f2name);
fprintf('    [VIZ] Saved: %s\n', f2name);

%% -----------------------------------------------------------------------
%  FIGURE 3: LS Channel Estimates across all pilot blocks
% -----------------------------------------------------------------------
fig3 = figure( ...
    'Name',     sprintf('LS Estimates | %s | Frame %d', cfg.CurrentMODCODName, frameIdx), ...
    'Position', [120 120 950 390], ...
    'Color',    'w');

pbIdx = (1:rxData.NumPilotBlocks)';

% Re(h_LS)
subplot(1,2,1);
plot(pbIdx, real(rxData.h_LS_perBlock), 'b-o', ...
    'LineWidth', 1.4, 'MarkerSize', 4, 'DisplayName', 'Re(h_{LS}) per block');
hold on;
yline(real(chData.h_tilde), 'r--', 'LineWidth', 2.0, ...
    'DisplayName', sprintf('True Re(h_{tilde}) = %.4f', real(chData.h_tilde)));
hold off;
grid on; box on;
xlabel('Pilot Block Index', 'FontSize', 10);
ylabel('Value',             'FontSize', 10);
title('Re(h_{LS})  vs  Pilot Block', 'FontSize', 11);
legend('Location', 'best', 'FontSize', 8);

% Im(h_LS)
subplot(1,2,2);
plot(pbIdx, imag(rxData.h_LS_perBlock), 'r-o', ...
    'LineWidth', 1.4, 'MarkerSize', 4, 'DisplayName', 'Im(h_{LS}) per block');
hold on;
yline(imag(chData.h_tilde), 'b--', 'LineWidth', 2.0, ...
    'DisplayName', sprintf('True Im(h_{tilde}) = %.4f', imag(chData.h_tilde)));
hold off;
grid on; box on;
xlabel('Pilot Block Index', 'FontSize', 10);
ylabel('Value',             'FontSize', 10);
title('Im(h_{LS})  vs  Pilot Block', 'FontSize', 11);
legend('Location', 'best', 'FontSize', 8);

sgtitle(sprintf( ...
    'LS Channel Estimates  |  MODCOD: %s  |  Frame %d  |  NMSE = %.4e', ...
    cfg.CurrentMODCODName, frameIdx, rxData.NMSE_LS), ...
    'FontSize', 12, 'FontWeight', 'bold');

f3name = fullfile(cfg.OutputDir, ...
    sprintf('fig3_ls_estimates_%s_fr%d.png', cfg.CurrentMODCODName, frameIdx));
saveas(fig3, f3name);
fprintf('    [VIZ] Saved: %s\n', f3name);

%% Summary
fprintf('    [VIZ] Frame %d  |  |h_tilde|=%.4f  |  A_p=%.2f dB  |  BER=%.5f  |  NMSE=%.4e\n', ...
    frameIdx, abs(chData.h_tilde), chData.rainAtten_dB, ...
    rxData.BER_LS, rxData.NMSE_LS);
end