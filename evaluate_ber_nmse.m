function evaluate_ber_nmse(cfg)
% =========================================================================
% EVALUATE_BER_NMSE — BER and NMSE vs SNR sweep plots
% =========================================================================
% Runs the full TX→CH→RX pipeline across cfg.SNR_test_range for each
% MODCOD and produces BER and NMSE plots comparable to Fig. 8 of the paper.
%
% Saves:
%   - BER_vs_SNR_LS.png
%   - NMSE_vs_SNR_LS.png
%   - results_<MODCOD>.csv  (numerical values)
% =========================================================================

snrRange = cfg.SNR_test_range;
numSNR   = length(snrRange);
numMOD   = length(cfg.MODCODs);

BER_LS_all   = zeros(numMOD, numSNR);
BER_noEq_all = zeros(numMOD, numSNR);
NMSE_LS_all  = zeros(numMOD, numSNR);

for modIdx = 1:numMOD
    cfg.CurrentMODCOD     = cfg.MODCODs(modIdx);
    cfg.CurrentMODCODName = cfg.MODCODNames{modIdx};

    fprintf('    [EVAL] MODCOD %s...\n', cfg.CurrentMODCODName);

    % Load ITU-R parameters once per MODCOD (same JSON reused)
    itu = load_itu_params(cfg.ParamsFile);

    % Generate one TX waveform (reused across all SNR points)
    txData = tx_dvbs2(cfg);

    for snrIdx = 1:numSNR
        snr_dB = snrRange(snrIdx);
        chData = channel_tropical(txData, cfg, snr_dB, itu);
        rxData = rx_dvbs2(txData, chData, cfg);

        BER_LS_all(modIdx, snrIdx)   = rxData.BER_LS;
        BER_noEq_all(modIdx, snrIdx) = rxData.BER_noEq;
        NMSE_LS_all(modIdx, snrIdx)  = rxData.NMSE_LS;
    end
end

%% -----------------------------------------------------------------------
%  PLOT 1: BER vs SNR
% -----------------------------------------------------------------------
colors  = {'b','r','g'};
markers = {'s','o','^'};

figure('Name','BER vs SNR — Tropical Channel (LS Estimator)', ...
    'Position',[100 100 800 560]);

for modIdx = 1:numMOD
    semilogy(snrRange, BER_LS_all(modIdx,:), ...
        [colors{modIdx} '-' markers{modIdx}], ...
        'LineWidth', 1.8, 'MarkerSize', 7, ...
        'DisplayName', sprintf('LS — %s', cfg.MODCODNames{modIdx}));
    hold on;
    semilogy(snrRange, BER_noEq_all(modIdx,:), ...
        [colors{modIdx} '--'], ...
        'LineWidth', 1.2, ...
        'DisplayName', sprintf('No Eq — %s', cfg.MODCODNames{modIdx}));
end

grid on;
xlabel('SNR (dB)', 'FontSize', 13);
ylabel('BER',      'FontSize', 13);
title({'BER vs SNR — Tropical Fading Channel (High Fading)', ...
    sprintf('DVB-S2 %s Frame | LS Channel Estimation | %.0f GHz Ka-band', ...
    cfg.FECFrame, cfg.Frequency_GHz)}, 'FontSize', 12);
legend('Location','southwest', 'FontSize', 10, 'NumColumns', 2);
ylim([1e-6 1]);
xlim([min(snrRange) max(snrRange)]);
saveas(gcf, fullfile(cfg.OutputDir, 'BER_vs_SNR_LS.png'));
fprintf('    [EVAL] Saved: BER_vs_SNR_LS.png\n');

%% -----------------------------------------------------------------------
%  PLOT 2: NMSE vs SNR
% -----------------------------------------------------------------------
figure('Name','NMSE vs SNR — LS Channel Estimator', ...
    'Position',[950 100 700 480]);

for modIdx = 1:numMOD
    semilogy(snrRange, NMSE_LS_all(modIdx,:), ...
        [colors{modIdx} '-' markers{modIdx}], ...
        'LineWidth', 1.8, 'MarkerSize', 7, ...
        'DisplayName', sprintf('LS — %s', cfg.MODCODNames{modIdx}));
    hold on;
end

grid on;
xlabel('SNR (dB)', 'FontSize', 13);
ylabel('NMSE',     'FontSize', 13);
title({'NMSE vs SNR — LS Channel Estimator', ...
    sprintf('Tropical Fading Channel | DVB-S2 %s Frame', cfg.FECFrame)}, ...
    'FontSize', 12);
legend('Location','northeast', 'FontSize', 10);
saveas(gcf, fullfile(cfg.OutputDir, 'NMSE_vs_SNR_LS.png'));
fprintf('    [EVAL] Saved: NMSE_vs_SNR_LS.png\n');

%% -----------------------------------------------------------------------
%  PRINT SUMMARY TABLE
% -----------------------------------------------------------------------
fprintf('\n========== BER/NMSE SUMMARY ==========\n');
for modIdx = 1:numMOD
    fprintf('\nMODCOD: %s\n', cfg.MODCODNames{modIdx});
    fprintf('%-10s %-15s %-15s %-15s\n', ...
        'SNR(dB)', 'BER_LS', 'BER_noEq', 'NMSE_LS');
    fprintf('%s\n', repmat('-',1,55));
    for snrIdx = 1:numSNR
        fprintf('%-10.1f %-15.6f %-15.6f %-15.6e\n', ...
            snrRange(snrIdx), ...
            BER_LS_all(modIdx,snrIdx), ...
            BER_noEq_all(modIdx,snrIdx), ...
            NMSE_LS_all(modIdx,snrIdx));
    end
end

%% -----------------------------------------------------------------------
%  SAVE NUMERICAL RESULTS TO CSV
% -----------------------------------------------------------------------
for modIdx = 1:numMOD
    name = cfg.MODCODNames{modIdx};
    T = table(snrRange(:), ...
              BER_LS_all(modIdx,:)', ...
              BER_noEq_all(modIdx,:)', ...
              NMSE_LS_all(modIdx,:)', ...
              'VariableNames', {'SNR_dB','BER_LS','BER_noEq','NMSE_LS'});
    fname = fullfile(cfg.OutputDir, sprintf('results_%s.csv', name));
    writetable(T, fname);
    fprintf('    [EVAL] Saved: results_%s.csv\n', name);
end

fprintf('\n    [EVAL] BER/NMSE evaluation complete.\n');
end
