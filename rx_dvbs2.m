function rxData = rx_dvbs2(txData, chData, cfg)

% RX_DVBS2 — DVB-S2 Receiver Module
% Performs matched filtering, pilot extraction, LS channel estimation,
% equalization, demodulation, and computes BER and NMSE.
%
% This module generates the INPUT VECTOR X_in for the BLSTM/GRU estimator
% as defined in Eq. (33) of Awad et al. (2023):
%
%   X_in = [ Re(y), Re(p), Re(h_LS), Im(y), Im(p), Im(h_LS) ]^T
%
% Where y = received pilot symbols, p = known pilot symbols, h_LS = LS estimate
% Each pilot block contributes ONE row to the dataset.
%
% INPUT:
%   txData  — struct from tx_dvbs2.m
%   chData  — struct from channel_tropical.m
%   cfg     — configuration struct
%
% OUTPUT:
%   rxData.X_in           — dataset matrix [NumPilotBlocks x 6] (Eq. 33)
%   rxData.h_LS_perBlock  — LS estimate per pilot block [NumPilotBlocks x 1]
%   rxData.h_true         — true channel coefficient (scalar)
%   rxData.BER_LS         — BER after LS equalization
%   rxData.BER_noEq       — BER without equalization
%   rxData.NMSE_LS        — NMSE of LS estimator
%   rxData.NumPilotBlocks — number of pilot blocks found
%   rxData.snr_dB         — SNR used
%   rxData.MODCOD         — MODCOD index
%   rxData.rainAtten_dB   — rain attenuation applied
% =========================================================================


% Matched Filter + Downsample to Symbol Rate
rxFilter = comm.RaisedCosineReceiveFilter( ...
    'RolloffFactor',         cfg.RolloffFactor, ...
    'FilterSpanInSymbols',   cfg.FilterSpanInSymbols, ...
    'InputSamplesPerSymbol', cfg.SamplesPerSymbol, ...
    'DecimationFactor',      cfg.SamplesPerSymbol);

rxSym_raw   = rxFilter(chData.rxWaveform);
filterDelay = cfg.FilterSpanInSymbols / 2;
rxSymbols   = rxSym_raw(filterDelay+1:end);
release(rxFilter);

% Align with TX symbols
txSymbols = txData.txSymbols;
minLen    = min(length(rxSymbols), length(txSymbols));
rxSymbols = rxSymbols(1:minLen);
txSymbols = txSymbols(1:minLen);

totalSymbols = length(rxSymbols);

% Pilot Block Detection and Extraction

% DVB-S2 Frame Structure (standard):
%   PLHEADER    : 90 symbols  (SOF + PLSCODE)
%   Data slot   : 90 symbols  (payload)
%   Pilot block : 36 symbols  (inserted every 16 data slots)
%
% Position of first pilot block:
%   PLHeader(90) + 16 slots * 90 symbols = 90 + 1440 = 1530
%
% Subsequent pilot blocks every: 36 (pilot) + 1440 (data) = 1476 symbols

pilotBlockLen       = cfg.PilotBlockLen;          % 36
slotsPerPilot       = cfg.SlotsPerPilot;          % 16
slotSize            = cfg.SlotSize;               % 90
PLHeaderLen         = cfg.PLHeaderLen;            % 90
symbolsPerDataChunk = slotsPerPilot * slotSize;   % 1440
pilotSpacing        = pilotBlockLen + symbolsPerDataChunk;  % 1476

% Known pilot value (DVB-S2 standard, un-modulated)
p_known = txData.pilotValue;  % (1+j)/sqrt(2)

% Find all pilot block start positions
pilotStartPos = [];
pos = PLHeaderLen + symbolsPerDataChunk + 1;  % 1-indexed MATLAB
while pos + pilotBlockLen - 1 <= totalSymbols
    pilotStartPos(end+1) = pos; %#ok<AGROW>
    pos = pos + pilotSpacing;
end

numPilotBlocks = length(pilotStartPos);

if numPilotBlocks == 0
    warning('[RX] No pilot blocks found. Check frame structure parameters.');
    rxData = build_empty_rxData(cfg, chData);
    return;
end

fprintf('    [RX] Pilot blocks found: %d\n', numPilotBlocks);

%  LS Channel Estimation Per Pilot Block (Eq. 12)

% For each pilot block b:
%   y_pilot = received pilot symbols  [36 x 1]
%   p_pilot = known pilot symbols     [36 x 1] (all equal to p_known)
%   h_LS_b  = mean(y_pilot ./ p_pilot)   scalar estimate for this block
%
% This gives us one h_LS estimate per pilot block.

% Preallocate
y_pilot_mean  = zeros(numPilotBlocks, 1);   % mean received pilot (complex)
p_pilot_mean  = p_known * ones(numPilotBlocks, 1);  % known (constant)
h_LS_perBlock = zeros(numPilotBlocks, 1);   % LS estimate per block

for pb = 1:numPilotBlocks
    startIdx = pilotStartPos(pb);
    endIdx   = startIdx + pilotBlockLen - 1;
    endIdx   = min(endIdx, totalSymbols);

    y_block  = rxSymbols(startIdx:endIdx);          % received pilots
    p_block  = p_known * ones(length(y_block), 1);  % known pilots

    % LS estimate (Eq. 12): h_LS = y / p (element-wise, then mean)
    h_LS_perBlock(pb) = mean(y_block ./ p_block);

    % Store mean received pilot for X_in construction
    y_pilot_mean(pb) = mean(y_block);
end

% Frame-level h_LS: average across all pilot blocks
h_LS_frame = mean(h_LS_perBlock);


% X_in = [ Re(y), Re(p), Re(h_LS), Im(y), Im(p), Im(h_LS) ]
%
% Dimensions: [numPilotBlocks x 6]
% Each row = one pilot block = one time step for BLSTM/GRU
%
% This is exactly what gets fed into the DL model in Python.

Re_y    = real(y_pilot_mean);
Re_p    = real(p_pilot_mean);
Re_hLS  = real(h_LS_perBlock);
Im_y    = imag(y_pilot_mean);
Im_p    = imag(p_pilot_mean);
Im_hLS  = imag(h_LS_perBlock);

X_in = [Re_y, Re_p, Re_hLS, Im_y, Im_p, Im_hLS];

% True channel label for each pilot block (target for DL training)
h_true = chData.h_true;
Re_h_true = real(h_true) * ones(numPilotBlocks, 1);
Im_h_true = imag(h_true) * ones(numPilotBlocks, 1);


% Equalize full received signal using frame-level h_LS
rxEq_LS    = rxSymbols / h_LS_frame;
rxEq_noEq  = rxSymbols;              % no equalization baseline

% Hard-decision QPSK demodulation
% QPSK Gray code: bit0 = sign(Re), bit1 = sign(Im)
demod_fn = @(s) double([real(s) >= 0, imag(s) >= 0]);

txBitPairs   = demod_fn(txSymbols);
rxBitPairs_LS   = demod_fn(rxEq_LS);
rxBitPairs_noEq = demod_fn(rxEq_noEq);

txB   = txBitPairs(:);
rxLS  = rxBitPairs_LS(:);
rxNoEq = rxBitPairs_noEq(:);

% Align lengths
minB = min([length(txB), length(rxLS), length(rxNoEq)]);
txB    = txB(1:minB);
rxLS   = rxLS(1:minB);
rxNoEq = rxNoEq(1:minB);


% BER and NMSE Computation
% BER
BER_LS   = sum(rxLS  ~= txB) / minB;
BER_noEq = sum(rxNoEq ~= txB) / minB;

% NMSE (Eq. 32 adapted):
%   NMSE = (1/K) * sum_k |h_true - h_LS_k|^2 / |h_true|^2
NMSE_LS = mean(abs(h_LS_perBlock - h_true).^2) / (abs(h_true)^2 + eps);

fprintf('    [RX] BER (LS)    = %.6f\n', BER_LS);
fprintf('    [RX] BER (no eq) = %.6f\n', BER_noEq);
fprintf('    [RX] NMSE (LS)   = %.6e\n', NMSE_LS);


rxData.X_in           = X_in;
rxData.h_LS_perBlock  = h_LS_perBlock;
rxData.h_LS_frame     = h_LS_frame;
rxData.h_true         = h_true;
rxData.Re_h_true      = Re_h_true;
rxData.Im_h_true      = Im_h_true;
rxData.BER_LS         = BER_LS;
rxData.BER_noEq       = BER_noEq;
rxData.NMSE_LS        = NMSE_LS;
rxData.NumPilotBlocks = numPilotBlocks;
rxData.snr_dB         = chData.snr_dB;
rxData.rainAtten_dB   = chData.rainAtten_dB;
rxData.MODCOD         = cfg.CurrentMODCOD;
rxData.MODCODName     = cfg.CurrentMODCODName;
rxData.pilotStartPos      = pilotStartPos;
rxData.PilotBlocksPerFrame = floor(numPilotBlocks / cfg.NumFrames);

end

%% Helper: empty struct on failure
function rxData = build_empty_rxData(cfg, chData)
    rxData.X_in           = [];
    rxData.h_LS_perBlock  = [];
    rxData.h_LS_frame     = NaN;
    rxData.h_true         = chData.h_true;
    rxData.Re_h_true      = [];
    rxData.Im_h_true      = [];
    rxData.BER_LS         = NaN;
    rxData.BER_noEq       = NaN;
    rxData.NMSE_LS        = NaN;
    rxData.NumPilotBlocks = 0;
    rxData.snr_dB         = chData.snr_dB;
    rxData.rainAtten_dB   = chData.rainAtten_dB;
    rxData.MODCOD         = cfg.CurrentMODCOD;
    rxData.MODCODName     = cfg.CurrentMODCODName;
    rxData.pilotStartPos  = [];
end