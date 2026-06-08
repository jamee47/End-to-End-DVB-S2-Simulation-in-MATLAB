function evaluate_ber_nmse(cfg, itu, allPilotPos)
% =========================================================================
% EVALUATE_BER_NMSE — BER and NMSE vs SNR for LS and MMSE estimators
% =========================================================================
% For each SNR point, generates cfg.NumEvalFrames independent frames,
% applies unique channel realisations, runs LS and MMSE estimation,
% and averages BER and NMSE over all frames.
%
% Produces plots comparable to Fig. 8 of Awad et al. (2023):
%   - NMSE vs SNR: LS, MMSE
%   - BER  vs SNR: LS, MMSE, No equalization, Theoretical QPSK
%
% INPUTS:
%   cfg          — configuration struct
%   itu          — struct from load_itu_params.m
%   allPilotPos  — struct of pilot positions per MODCOD
% =========================================================================

snrRange      = cfg.SNR_test_range;
numSNR        = length(snrRange);
numMOD        = length(cfg.MODCODs);
numEvalFrames = cfg.NumEvalFrames;

fprintf('\n>>> BER/NMSE evaluation: %d MODCODs x %d SNR pts x %d frames/pt\n', ...
    numMOD, numSNR, numEvalFrames);

% Preallocate result arrays
BER_LS_all   = zeros(numMOD, numSNR);
BER_MMSE_all = zeros(numMOD, numSNR);
BER_noEq_all = zeros(numMOD, numSNR);
NMSE_LS_all  = zeros(numMOD, numSNR);
NMSE_MMSE_all= zeros(numMOD, numSNR);

for modIdx = 1:numMOD
    cfg.CurrentMODCOD     = cfg.MODCODs(modIdx);
    cfg.CurrentMODCODName = cfg.MODCODNames{modIdx};
    pilotPos              = allPilotPos.(cfg.CurrentMODCODName);

    fprintf('\n    [EVAL] MODCOD: %s\n', cfg.CurrentMODCODName);
    fprintf('    %-8s %-12s %-12s %-12s %-12s %-12s\n', ...
        'SNR(dB)','BER_LS','BER_MMSE','BER_noEq','NMSE_LS','NMSE_MMSE');
    fprintf('    %s\n', repmat('-',1,68));

    for snrIdx = 1:numSNR
        snr_dB = snrRange(snrIdx);

        % Accumulate over frames
        berLS_sum    = 0;
        berMMSE_sum  = 0;
        berNoEq_sum  = 0;
        nmseLS_sum   = 0;
        nmseMMSE_sum = 0;
        validFrames  = 0;

        for f = 1:numEvalFrames
            try
                txFrame = tx_dvbs2_frame(cfg);
                chFrame = channel_tropical(txFrame, cfg, snr_dB, itu);
                rxFrame = rx_dvbs2_frame(txFrame, chFrame, cfg, pilotPos);

                % Skip frames where BER is NaN (FEC decode failed)
                if isnan(rxFrame.BER_LS), continue; end

                berLS_sum    = berLS_sum    + rxFrame.BER_LS;
                berMMSE_sum  = berMMSE_sum  + rxFrame.BER_MMSE;
                berNoEq_sum  = berNoEq_sum  + rxFrame.BER_noEq;
                nmseLS_sum   = nmseLS_sum   + rxFrame.NMSE_LS;
                nmseMMSE_sum = nmseMMSE_sum + rxFrame.NMSE_MMSE;
                validFrames  = validFrames  + 1;
            catch ME
                % Skip frames that error (e.g. FEC sync failure)
                continue;
            end
        end

        if validFrames == 0
            BER_LS_all(modIdx,snrIdx)    = NaN;
            BER_MMSE_all(modIdx,snrIdx)  = NaN;
            BER_noEq_all(modIdx,snrIdx)  = NaN;
            NMSE_LS_all(modIdx,snrIdx)   = NaN;
            NMSE_MMSE_all(modIdx,snrIdx) = NaN;
        else
            BER_LS_all(modIdx,snrIdx)    = berLS_sum    / validFrames;
            BER_MMSE_all(modIdx,snrIdx)  = berMMSE_sum  / validFrames;
            BER_noEq_all(modIdx,snrIdx)  = berNoEq_sum  / validFrames;
            NMSE_LS_all(modIdx,snrIdx)   = nmseLS_sum   / validFrames;
            NMSE_MMSE_all(modIdx,snrIdx) = nmseMMSE_sum / validFrames;
        end

        fprintf('    %-8.1f %-12.6f %-12.6f %-12.6f %-12.4e %-12.4e\n', ...
            snr_dB, ...
            BER_LS_all(modIdx,snrIdx), ...
            BER_MMSE_all(modIdx,snrIdx), ...
            BER_noEq_all(modIdx,snrIdx), ...
            NMSE_LS_all(modIdx,snrIdx), ...
            NMSE_MMSE_all(modIdx,snrIdx));
    end
end

%% -----------------------------------------------------------------------
%  PLOT 1: NMSE vs SNR
%  Comparable to Fig. 8(a)(c)(e) of Awad et al. 2023
% -----------------------------------------------------------------------
colors  = {'b','r','g','m'};
markers = {'s','o','^','d'};
lw = 1.8; ms = 7;

fig1 = figure('Name','NMSE vs SNR', 'Position',[100 200 820 520], 'Color','w');

for modIdx = 1:numMOD
    c = colors{modIdx}; mk = markers{modIdx};
    name = cfg.MODCODNames{modIdx};

    semilogy(snrRange, NMSE_LS_all(modIdx,:), ...
        [c '--' mk], 'LineWidth',lw, 'MarkerSize',ms, ...
        'DisplayName', sprintf('LS — %s', name));
    hold on;
    semilogy(snrRange, NMSE_MMSE_all(modIdx,:), ...
        [c '-' mk], 'LineWidth',lw, 'MarkerSize',ms, ...
        'DisplayName', sprintf('MMSE — %s (upper bound)', name));
end

grid on;
xlabel('E_s/N_0  (dB)', 'FontSize',13);
ylabel('NMSE',          'FontSize',13);
title({'NMSE vs SNR — Tropical Fading Channel', ...
    sprintf('DVB-S2 %s Frame | %.0f GHz | %d frames/point', ...
    cfg.FECFrame, cfg.Frequency_GHz, numEvalFrames)}, 'FontSize',12);
legend('Location','northeast','FontSize',10);
xlim([min(snrRange) max(snrRange)]);

p1 = fullfile(cfg.OutputDir,'NMSE_vs_SNR.png');
saveas(fig1, p1);
fprintf('\n    [EVAL] Saved: %s\n', p1);

%% -----------------------------------------------------------------------
%  PLOT 2: BER vs SNR
%  Comparable to Fig. 8(b)(d)(f) of Awad et al. 2023
% -----------------------------------------------------------------------
fig2 = figure('Name','BER vs SNR', 'Position',[130 130 820 520], 'Color','w');

% Theoretical QPSK BER in AWGN (reference line)
snr_lin  = 10.^(snrRange/10);
BER_theo = 0.5 * erfc(sqrt(snr_lin));

for modIdx = 1:numMOD
    c = colors{modIdx}; mk = markers{modIdx};
    name = cfg.MODCODNames{modIdx};

    % No equalization
    semilogy(snrRange, BER_noEq_all(modIdx,:), ...
        [c ':'], 'LineWidth',1.2, ...
        'DisplayName', sprintf('No Eq — %s', name));
    hold on;

    % LS
    semilogy(snrRange, BER_LS_all(modIdx,:), ...
        [c '--' mk], 'LineWidth',lw, 'MarkerSize',ms, ...
        'DisplayName', sprintf('LS — %s', name));

    % MMSE (upper bound)
    semilogy(snrRange, BER_MMSE_all(modIdx,:), ...
        [c '-' mk], 'LineWidth',lw, 'MarkerSize',ms, ...
        'DisplayName', sprintf('MMSE — %s', name));
end

% Perfect CSI reference
semilogy(snrRange, BER_theo, 'k-.', 'LineWidth',1.5, ...
    'DisplayName','Theoretical QPSK (AWGN)');

grid on;
xlabel('E_s/N_0  (dB)', 'FontSize',13);
ylabel('BER',            'FontSize',13);
title({'BER vs SNR — Tropical Fading Channel', ...
    sprintf('DVB-S2 %s Frame | %.0f GHz | %d frames/point', ...
    cfg.FECFrame, cfg.Frequency_GHz, numEvalFrames)}, 'FontSize',12);
legend('Location','southwest','FontSize',10,'NumColumns',2);
ylim([1e-5 1]);
xlim([min(snrRange) max(snrRange)]);

p2 = fullfile(cfg.OutputDir,'BER_vs_SNR.png');
saveas(fig2, p2);
fprintf('    [EVAL] Saved: %s\n', p2);

%% -----------------------------------------------------------------------
%  Save numerical results to CSV
% -----------------------------------------------------------------------
for modIdx = 1:numMOD
    name = cfg.MODCODNames{modIdx};
    T = table(snrRange(:), ...
        BER_LS_all(modIdx,:)',    ...
        BER_MMSE_all(modIdx,:)',  ...
        BER_noEq_all(modIdx,:)',  ...
        NMSE_LS_all(modIdx,:)',   ...
        NMSE_MMSE_all(modIdx,:)', ...
        'VariableNames', ...
        {'SNR_dB','BER_LS','BER_MMSE','BER_noEq','NMSE_LS','NMSE_MMSE'});

    fname = fullfile(cfg.OutputDir, sprintf('results_%s.csv', name));
    writetable(T, fname);
    fprintf('    [EVAL] Saved: results_%s.csv\n', name);
end

fprintf('\n    [EVAL] Complete.\n');
end
