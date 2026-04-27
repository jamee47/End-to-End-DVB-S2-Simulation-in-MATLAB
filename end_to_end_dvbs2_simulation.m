clc, clearvars, close all;
% Author : Mohtasim Al Jamee
% Modified: Added LS Channel Estimation and BER vs SNR sweep

%% DVB-S2 Configuration
wavegen = dvbs2WaveformGenerator;
wavegen.StreamFormat     = "TS";
wavegen.NumInputStreams  = 1;
wavegen.FECFrame         = "normal";
wavegen.MODCOD           = 10;          % QPSK 1/2
wavegen.DFL              = 15928;
wavegen.HasPilots        = 1;          % Pilots MUST be on for LS estimation
wavegen.RolloffFactor    = 0.35;
wavegen.FilterSpanInSymbols = 10;
wavegen.SamplesPerSymbol = 4;

%% Input Bit Source
NumPLFrames = 10;   % Increased for reliable BER statistics
pn = comm.PNSequence('Polynomial','x9+x5+1', ...
    'InitialConditions',[zeros(1,8) 1], ...
    'VariableSizeOutput',true, ...
    'MaximumOutputSize',[64800*NumPLFrames 1]);

% Build transport stream packets
syncBits = [0;1;0;0;0;1;1;1];
pktLen   = 1496;
data     = cell(1, wavegen.NumInputStreams);
for i = 1:wavegen.NumInputStreams
    numPkts  = wavegen.MinNumPackets(i) * NumPLFrames;
    reset(pn);
    rawpkts  = pn(pktLen * numPkts);
    txRawPkts = reshape(rawpkts, pktLen, numPkts);
    txPkts   = [repmat(syncBits,1,numPkts); txRawPkts];
    data{i}  = txPkts(:);
end

% Store original transmitted bits for BER calculation
txBits = data{1};

%% Generate Clean Waveform (no noise yet — we add noise per SNR point)
waveform_clean = [wavegen(data); flushFilter(wavegen)];
Fs = 1e6;  % Sample rate in Hz

%% -------------------------------------------------------------------
%  LS CHANNEL ESTIMATION + BER vs SNR SWEEP
% -------------------------------------------------------------------
% Explanation of LS estimation:
%   The received signal model is:  y = h * x + n
%   At pilot positions, x = p (known), so:
%   h_LS = y_pilot / p_pilot   (element-wise division)
%   The estimated h_LS is then used to equalize:
%   x_hat = y / h_LS
%
% For AWGN: h = 1 (no fading), so h_LS should estimate ~1.
% This gives us a baseline before adding fading later.
% -------------------------------------------------------------------

SNR_dB_range = 0:2:30;           % SNR sweep from 0 to 30 dB
numSNR       = length(SNR_dB_range);

BER_no_estimation   = zeros(1, numSNR);   % Ideal: no channel correction
BER_LS_estimation   = zeros(1, numSNR);   % With LS channel estimation
BER_perfect_CSI     = zeros(1, numSNR);   % Perfect channel knowledge (upper bound)

% DVB-S2 demodulator configuration (must match transmitter)
demod = comm.OFDMDemodulator;   % placeholder — we use dvbs2 receiver below

% Use MATLAB's dvbs2 receiver-side objects
rxConfig.MODCOD         = wavegen.MODCOD;
rxConfig.FECFrame       = wavegen.FECFrame;
rxConfig.HasPilots      = wavegen.HasPilots;
rxConfig.RolloffFactor  = wavegen.RolloffFactor;
rxConfig.SamplesPerSymbol = wavegen.SamplesPerSymbol;
rxConfig.FilterSpanInSymbols = wavegen.FilterSpanInSymbols;

fprintf('Running BER vs SNR sweep...\n');
fprintf('%-10s %-20s %-20s %-20s\n','SNR(dB)','BER(No Est.)','BER(LS Est.)','BER(Perfect CSI)');
fprintf('%s\n', repmat('-',1,72));

for snrIdx = 1:numSNR
    snr_dB = SNR_dB_range(snrIdx);

    %% Step 1: Add AWGN to clean waveform
    rxWaveform = awgn(waveform_clean, snr_dB, 'measured');

    %% Step 2: Extract pilot symbols for LS estimation
    % In DVB-S2, pilots are known: I = 1/sqrt(2), Q = 1/sqrt(2)
    % Pilot blocks are inserted every 16 slots (each slot = 90 symbols)
    % Each pilot block = 36 symbols
    % We work at symbol level (after matched filtering / downsampling)

    % --- Matched filter + downsample to 1 sample/symbol ---
    % Design receive square-root raised cosine filter (matches transmitter)
    rxFilter = comm.RaisedCosineReceiveFilter( ...
        'RolloffFactor',        wavegen.RolloffFactor, ...
        'FilterSpanInSymbols',  wavegen.FilterSpanInSymbols, ...
        'InputSamplesPerSymbol',wavegen.SamplesPerSymbol, ...
        'DecimationFactor',     wavegen.SamplesPerSymbol);

    rxSymbols_raw = rxFilter(rxWaveform);
    release(rxFilter);

    % Remove filter delay (half filter span)
    filterDelay = wavegen.FilterSpanInSymbols / 2;
    rxSymbols   = rxSymbols_raw(filterDelay+1:end);

    % --- Regenerate clean (noiseless) symbols for pilot reference ---
    txFilter = comm.RaisedCosineReceiveFilter( ...
        'RolloffFactor',        wavegen.RolloffFactor, ...
        'FilterSpanInSymbols',  wavegen.FilterSpanInSymbols, ...
        'InputSamplesPerSymbol',wavegen.SamplesPerSymbol, ...
        'DecimationFactor',     wavegen.SamplesPerSymbol);

    txSymbols_raw = txFilter(waveform_clean);
    release(txFilter);
    txSymbols = txSymbols_raw(filterDelay+1:end);

    % Align lengths
    minLen    = min(length(rxSymbols), length(txSymbols));
    rxSymbols = rxSymbols(1:minLen);
    txSymbols = txSymbols(1:minLen);

    %% Step 3: LS Channel Estimation
    % DVB-S2 pilot symbol value (known constant)
    pilotValue = (1/sqrt(2)) + 1j*(1/sqrt(2));

    % Pilot block structure:
    %   - PLHEADER: 90 symbols (slots 1)
    %   - Data slots: 16 slots × 90 symbols = 1440 symbols between pilots
    %   - Pilot block: 36 symbols
    % Total symbols per PLFRAME (QPSK normal):
    %   PLHEADER(90) + 360 data slots × 90 + pilot blocks × 36
    % For MODCOD=1 (QPSK 1/2 normal): 33282 symbols per PLFRAME approx.

    slotSize      = 90;
    pilotBlockLen = 36;
    slotsPerPilot = 16;
    symbolsPerDataChunk = slotsPerPilot * slotSize;  % 1440 symbols between pilots

    % Extract pilot positions and estimate channel
    totalSymbols  = length(rxSymbols);
    pilotPositions = [];
    pos = 90 + symbolsPerDataChunk;  % skip PLHEADER, first data chunk
    while pos + pilotBlockLen - 1 <= totalSymbols
        pilotPositions(end+1) = pos; %#ok<AGROW>
        pos = pos + pilotBlockLen + symbolsPerDataChunk;
    end

    numPilotBlocks = length(pilotPositions);

    if numPilotBlocks == 0
        % Fallback: no pilots found, h_LS = 1
        h_LS = 1;
    else
        h_estimates = zeros(numPilotBlocks, 1);
        for pb = 1:numPilotBlocks
            idx = pilotPositions(pb) : pilotPositions(pb) + pilotBlockLen - 1;
            idx = idx(idx <= totalSymbols);
            if isempty(idx), continue; end
            % LS estimate: h_LS = mean(y_pilot / p_pilot)
            h_estimates(pb) = mean(rxSymbols(idx) / pilotValue);
        end
        % Average h_LS over all pilot blocks in the frame
        h_LS = mean(h_estimates);
    end

    %% Step 4: Equalization using h_LS
    % Equalized signal: x_hat = y / h_LS
    rxEqualized_LS = rxSymbols / h_LS;

    % Perfect CSI equalization (h=1 for AWGN, but computed from tx symbols)
    h_perfect      = mean(txSymbols(txSymbols ~= 0)) / ...
                     mean(rxSymbols(txSymbols ~= 0) ./ txSymbols(txSymbols ~= 0));
    % For pure AWGN h_perfect = 1; keep it explicit for when fading is added
    rxEqualized_perfect = rxSymbols / 1;  % h=1 for AWGN

    %% Step 5: QPSK Demodulation (hard decision)
    % QPSK Gray-coded: decision based on sign of I and Q
    demodulate = @(sig) [real(sig) >= 0, imag(sig) >= 0];

    % Reference: demodulate clean tx symbols to get tx bit pairs
    txBitPairs_ref = demodulate(txSymbols);
    txBits_sym     = txBitPairs_ref(:);

    % No estimation (raw received, no equalization)
    rxBits_noEst   = demodulate(rxSymbols);
    rxBits_noEst   = rxBits_noEst(:);

    % LS estimation
    rxBits_LS      = demodulate(rxEqualized_LS);
    rxBits_LS      = rxBits_LS(:);

    % Perfect CSI
    rxBits_perfect = demodulate(rxEqualized_perfect);
    rxBits_perfect = rxBits_perfect(:);

    % Align lengths for BER
    minBits = min([length(txBits_sym), length(rxBits_noEst), ...
                   length(rxBits_LS),  length(rxBits_perfect)]);

    txB = txBits_sym(1:minBits);

    %% Step 6: Compute BER
    BER_no_estimation(snrIdx) = sum(rxBits_noEst(1:minBits)  ~= txB) / minBits;
    BER_LS_estimation(snrIdx) = sum(rxBits_LS(1:minBits)     ~= txB) / minBits;
    BER_perfect_CSI(snrIdx)   = sum(rxBits_perfect(1:minBits) ~= txB) / minBits;

    fprintf('%-10.1f %-20.6f %-20.6f %-20.6f\n', snr_dB, ...
        BER_no_estimation(snrIdx), BER_LS_estimation(snrIdx), BER_perfect_CSI(snrIdx));
end

%% -------------------------------------------------------------------
%  PLOTS
% -------------------------------------------------------------------

figure('Name','BER vs SNR - DVB-S2 with LS Channel Estimation', ...
       'Position',[100 100 800 550]);

semilogy(SNR_dB_range, BER_no_estimation,  'r--o', 'LineWidth', 1.8, ...
    'MarkerSize', 7, 'DisplayName', 'No Equalization');
hold on;
semilogy(SNR_dB_range, BER_LS_estimation,  'b-s',  'LineWidth', 1.8, ...
    'MarkerSize', 7, 'DisplayName', 'LS Channel Estimation');
semilogy(SNR_dB_range, BER_perfect_CSI,    'k-^',  'LineWidth', 1.8, ...
    'MarkerSize', 7, 'DisplayName', 'Perfect CSI (upper bound)');

grid on;
xlabel('SNR (dB)',    'FontSize', 13);
ylabel('BER',         'FontSize', 13);
title({'DVB-S2 BER vs SNR', ...
       sprintf('MODCOD %d | QPSK 1/2 | Normal Frame | AWGN Channel', wavegen.MODCOD)}, ...
       'FontSize', 13);
legend('Location','southwest', 'FontSize', 11);
ylim([1e-6 1]);
xlim([0 30]);

% Add theoretical QPSK BER curve for reference
snr_lin  = 10.^(SNR_dB_range/10);
BER_theo = qfunc(sqrt(2 * snr_lin));
semilogy(SNR_dB_range, BER_theo, 'm:', 'LineWidth', 1.5, ...
    'DisplayName', 'Theoretical QPSK (AWGN)');
legend('Location','southwest', 'FontSize', 11);

hold off;

%% -------------------------------------------------------------------
%  NMSE of LS Estimator
% -------------------------------------------------------------------
% For AWGN, true h = 1. NMSE = |h_LS - h_true|^2 / |h_true|^2
% We recompute per SNR for the NMSE plot

fprintf('\nComputing NMSE of LS estimator...\n');
NMSE_LS = zeros(1, numSNR);
h_true  = 1;  % AWGN: no fading, true channel = 1

for snrIdx = 1:numSNR
    snr_dB     = SNR_dB_range(snrIdx);
    rxWaveform = awgn(waveform_clean, snr_dB, 'measured');

    rxFilt = comm.RaisedCosineReceiveFilter( ...
        'RolloffFactor',        wavegen.RolloffFactor, ...
        'FilterSpanInSymbols',  wavegen.FilterSpanInSymbols, ...
        'InputSamplesPerSymbol',wavegen.SamplesPerSymbol, ...
        'DecimationFactor',     wavegen.SamplesPerSymbol);
    rxSym_raw  = rxFilt(rxWaveform);
    release(rxFilt);
    filterDelay = wavegen.FilterSpanInSymbols / 2;
    rxSym       = rxSym_raw(filterDelay+1:end);
    totalSym    = length(rxSym);

    pilotPos = [];
    pos = 90 + symbolsPerDataChunk;
    while pos + pilotBlockLen - 1 <= totalSym
        pilotPos(end+1) = pos; %#ok<AGROW>
        pos = pos + pilotBlockLen + symbolsPerDataChunk;
    end

    if isempty(pilotPos)
        NMSE_LS(snrIdx) = NaN;
        continue;
    end

    h_est_all = zeros(length(pilotPos),1);
    for pb = 1:length(pilotPos)
        idx = pilotPos(pb) : pilotPos(pb) + pilotBlockLen - 1;
        idx = idx(idx <= totalSym);
        if isempty(idx), continue; end
        h_est_all(pb) = mean(rxSym(idx) / pilotValue);
    end
    h_LS_mean = mean(h_est_all);
    NMSE_LS(snrIdx) = abs(h_LS_mean - h_true)^2 / abs(h_true)^2;
end

figure('Name','NMSE of LS Channel Estimator vs SNR', ...
       'Position',[950 100 700 480]);
semilogy(SNR_dB_range, NMSE_LS, 'b-s', 'LineWidth', 1.8, ...
    'MarkerSize', 7, 'DisplayName', 'LS Estimator NMSE');
grid on;
xlabel('SNR (dB)', 'FontSize', 13);
ylabel('NMSE',     'FontSize', 13);
title('NMSE of LS Channel Estimator vs SNR (AWGN)', 'FontSize', 13);
legend('Location','northeast', 'FontSize', 11);

%% -------------------------------------------------------------------
%  Summary Table
% -------------------------------------------------------------------
fprintf('\n========== SUMMARY TABLE ==========\n');
fprintf('%-10s %-20s %-20s %-20s %-15s\n', ...
    'SNR(dB)','BER(No Est.)','BER(LS Est.)','BER(Perfect)','NMSE(LS)');
fprintf('%s\n', repmat('-',1,87));
for i = 1:numSNR
    fprintf('%-10.1f %-20.6f %-20.6f %-20.6f %-15.6e\n', ...
        SNR_dB_range(i), BER_no_estimation(i), ...
        BER_LS_estimation(i), BER_perfect_CSI(i), NMSE_LS(i));
end