function rxData = rx_dvbs2_frame(txData, chData, cfg, pilotPos)
% =========================================================================
% RX_DVBS2_FRAME — Receiver for ONE DVB-S2 frame
% =========================================================================
% Implements:
%   - Matched filtering (SRRC)
%   - LS channel estimation per pilot block (Eq. 12)
%   - MMSE channel estimation per pilot block (Eq. 17)
%   - Equalization by both estimators
%   - Pre-FEC BER via QPSK symbol error rate on equalized symbols
%   - NMSE per Eq. 32
%   - X_in feature matrix per Eq. 33
%
% NOTE on BER:
%   MATLAB has no standalone dvbs2Receiver object. The paper's BER
%   is the pre-FEC symbol error rate after equalization, mapped to
%   bit errors assuming Gray-coded QPSK (2 bits/symbol).
%   This is consistent with how simulation-based papers measure BER
%   before the FEC coding gain is applied.
%   BER = symbol_errors / total_symbols  (Gray-coded QPSK: 1 error/sym)
%
% NOTE on MMSE:
%   For flat slow-fading (constant h within one frame):
%   h_MMSE = W * h_LS,  W = R_hh / (R_hh + sigma_n2/Np)
%   where R_hh = mean channel power, sigma_n2 = noise power per symbol,
%   Np = 36 (pilot block averaging gain).
%
% INPUT:
%   txData   — struct from tx_dvbs2_frame.m
%   chData   — struct from channel_tropical.m
%   cfg      — configuration struct
%   pilotPos — pilot positions from detect_pilot_positions.m
%
% OUTPUT:
%   rxData.X_in            — [numPilotBlocks x 6]  feature matrix
%   rxData.h_LS_perBlock   — complex LS estimate per pilot block
%   rxData.h_MMSE_perBlock — complex MMSE estimate per pilot block
%   rxData.h_LS_frame      — frame-level LS  (mean over blocks)
%   rxData.h_MMSE_frame    — frame-level MMSE (mean over blocks)
%   rxData.BER_LS          — pre-FEC BER after LS equalization
%   rxData.BER_MMSE        — pre-FEC BER after MMSE equalization
%   rxData.BER_noEq        — pre-FEC BER without equalization
%   rxData.NMSE_LS         — NMSE of LS estimator   (Eq. 32)
%   rxData.NMSE_MMSE       — NMSE of MMSE estimator (Eq. 32)
%   rxData.NumPilotBlocks  — number of pilot blocks found
% =========================================================================

pilotBlockLen = cfg.PilotBlockLen;   % 36 symbols per pilot block
p_known       = txData.pilotValue;   % (1+j)/sqrt(2)

%% -----------------------------------------------------------------------
%  STEP 1: Matched filter — downsample to 1 sample/symbol
% -----------------------------------------------------------------------
rxFilter = comm.RaisedCosineReceiveFilter( ...
    'RolloffFactor',         cfg.RolloffFactor, ...
    'FilterSpanInSymbols',   cfg.FilterSpanInSymbols, ...
    'InputSamplesPerSymbol', cfg.SamplesPerSymbol, ...
    'DecimationFactor',      cfg.SamplesPerSymbol);

rxSym_raw   = rxFilter(chData.rxWaveform);
filterDelay = cfg.FilterSpanInSymbols / 2;
rxSymbols   = rxSym_raw(filterDelay+1:end);
release(rxFilter);

txSymbols = txData.txSymbols;

% Align lengths
minLen    = min(length(rxSymbols), length(txSymbols));
rxSymbols = rxSymbols(1:minLen);
txSymbols = txSymbols(1:minLen);
totalSyms = minLen;

numPilotBlocks = length(pilotPos);

%% -----------------------------------------------------------------------
%  STEP 2: LS estimation per pilot block (Eq. 12)
%
%  h_LS_b = mean(y_b / p)
%  Averaging Np=36 symbols reduces noise variance by factor Np.
% -----------------------------------------------------------------------
y_pilot_mean  = zeros(numPilotBlocks, 1);
p_pilot_mean  = p_known * ones(numPilotBlocks, 1);
h_LS_perBlock = zeros(numPilotBlocks, 1);

for pb = 1:numPilotBlocks
    startIdx = pilotPos(pb);
    endIdx   = min(startIdx + pilotBlockLen - 1, totalSyms);

    if startIdx > totalSyms
        % Pad with last known estimate if frame is shorter than expected
        h_LS_perBlock(pb) = h_LS_perBlock(max(pb-1,1));
        y_pilot_mean(pb)  = y_pilot_mean(max(pb-1,1));
        continue;
    end

    y_block           = rxSymbols(startIdx:endIdx);
    h_LS_perBlock(pb) = mean(y_block ./ p_known);
    y_pilot_mean(pb)  = mean(y_block);
end

h_LS_frame = mean(h_LS_perBlock);

%% -----------------------------------------------------------------------
%  STEP 3: MMSE estimation per pilot block (Eq. 17, simplified scalar)
%
%  For flat slow-fading:
%    h_MMSE = W * h_LS
%    W = R_hh / (R_hh + sigma_n2 / Np)
%
%  R_hh  = E[|h|^2] estimated from pilot blocks (bias-corrected)
%  sigma_n2 = noise power per symbol = |p|^2 / SNR_linear
%  Np = pilotBlockLen = 36 (averaging gain already in h_LS_perBlock)
% -----------------------------------------------------------------------
snr_lin  = 10^(chData.snr_dB / 10);
p_power  = abs(p_known)^2;            % = 0.5 for QPSK pilot
sigma_n2 = p_power / snr_lin;         % noise variance per symbol

% Estimate channel power: R_hh = E[|h_LS|^2] - sigma_n2/Np
Np   = pilotBlockLen;
R_hh = mean(abs(h_LS_perBlock).^2) - sigma_n2 / Np;
R_hh = max(real(R_hh), 1e-12);        % must be positive

% Wiener filter coefficient
W_mmse = R_hh / (R_hh + sigma_n2 / Np);

% Apply scalar MMSE filter to all pilot block estimates
h_MMSE_perBlock = W_mmse * h_LS_perBlock;
h_MMSE_frame    = mean(h_MMSE_perBlock);

%% -----------------------------------------------------------------------
%  STEP 4: X_in feature matrix (Eq. 33)
%  [Re(y), Re(p), Re(h_LS), Im(y), Im(p), Im(h_LS)] per pilot block
% -----------------------------------------------------------------------
X_in = [real(y_pilot_mean), real(p_pilot_mean), real(h_LS_perBlock), ...
        imag(y_pilot_mean), imag(p_pilot_mean), imag(h_LS_perBlock)];

%% -----------------------------------------------------------------------
%  STEP 5: Equalization
%  y_eq = y / h_est  (scalar division — flat fading assumption)
% -----------------------------------------------------------------------
rxEq_LS   = rxSymbols / h_LS_frame;
rxEq_MMSE = rxSymbols / h_MMSE_frame;
rxEq_noEq = rxSymbols;                 % no correction, h_est = 1

%% -----------------------------------------------------------------------
%  STEP 6: Pre-FEC BER via QPSK Gray-coded symbol decisions
%
%  QPSK decision regions (Gray code):
%    bit0 = (Re(s) >= 0) → maps to I component
%    bit1 = (Im(s) >= 0) → maps to Q component
%
%  For data symbols only (exclude pilot positions from BER calculation)
%  Pilots are known and always correct — including them inflates BER.
% -----------------------------------------------------------------------

% Build data symbol mask — exclude pilot positions
dataMask = true(totalSyms, 1);
for pb = 1:numPilotBlocks
    startIdx = pilotPos(pb);
    endIdx   = min(startIdx + pilotBlockLen - 1, totalSyms);
    dataMask(startIdx:endIdx) = false;
end
% Also exclude PLHEADER
plhLen = cfg.PLHeaderLen;
if plhLen < totalSyms
    dataMask(1:plhLen) = false;
end

% Apply mask
txData_syms   = txSymbols(dataMask);
rxData_LS     = rxEq_LS(dataMask);
rxData_MMSE   = rxEq_MMSE(dataMask);
rxData_noEq   = rxEq_noEq(dataMask);

if isempty(txData_syms)
    % Fallback: use all symbols if mask removes everything
    txData_syms = txSymbols;
    rxData_LS   = rxEq_LS;
    rxData_MMSE = rxEq_MMSE;
    rxData_noEq = rxEq_noEq;
end

% Gray-coded QPSK demodulation
% TX symbols are the reference — compare sign of I and Q components
qpsk_decide = @(s) [real(s) >= 0, imag(s) >= 0];

txBits_sym  = qpsk_decide(txData_syms);   txBits_sym  = txBits_sym(:);
rxBits_LS   = qpsk_decide(rxData_LS);     rxBits_LS   = rxBits_LS(:);
rxBits_MMSE = qpsk_decide(rxData_MMSE);   rxBits_MMSE = rxBits_MMSE(:);
rxBits_noEq = qpsk_decide(rxData_noEq);   rxBits_noEq = rxBits_noEq(:);

% Align lengths
nb = min([length(txBits_sym), length(rxBits_LS), ...
          length(rxBits_MMSE), length(rxBits_noEq)]);

BER_LS   = sum(txBits_sym(1:nb) ~= rxBits_LS(1:nb))   / nb;
BER_MMSE = sum(txBits_sym(1:nb) ~= rxBits_MMSE(1:nb)) / nb;
BER_noEq = sum(txBits_sym(1:nb) ~= rxBits_noEq(1:nb)) / nb;

%% -----------------------------------------------------------------------
%  STEP 7: NMSE (Eq. 32)
%
%  NMSE = mean_b(|h_true - h_est_b|^2) / |h_true|^2
%
%  h_true = h_tilde (scalar for this frame — same channel for all blocks)
%  Averaging over pilot blocks gives a per-frame NMSE estimate.
% -----------------------------------------------------------------------
h_true   = chData.h_tilde;
h_power  = abs(h_true)^2 + eps;

NMSE_LS   = mean(abs(h_LS_perBlock   - h_true).^2) / h_power;
NMSE_MMSE = mean(abs(h_MMSE_perBlock - h_true).^2) / h_power;

%% Pack output
rxData.X_in            = X_in;
rxData.h_LS_perBlock   = h_LS_perBlock;
rxData.h_MMSE_perBlock = h_MMSE_perBlock;
rxData.h_LS_frame      = h_LS_frame;
rxData.h_MMSE_frame    = h_MMSE_frame;
rxData.BER_LS          = BER_LS;
rxData.BER_MMSE        = BER_MMSE;
rxData.BER_noEq        = BER_noEq;
rxData.NMSE_LS         = NMSE_LS;
rxData.NMSE_MMSE       = NMSE_MMSE;
rxData.NumPilotBlocks  = numPilotBlocks;
rxData.rxSymbols       = rxSymbols;
rxData.txSymbols       = txSymbols;
rxData.W_mmse          = W_mmse;
rxData.R_hh            = R_hh;
rxData.sigma_n2        = sigma_n2;
end