function txData = tx_dvbs2_frame(cfg)
% =========================================================================
% TX_DVB2_FRAME — Generate exactly ONE DVB-S2 physical layer frame
% =========================================================================
% Called once per frame in the dataset generation loop.
% Each call generates unique random bits → unique encoded waveform.
% The channel is applied separately in channel_tropical.m.
%
% OUTPUT:
%   txData.waveform        — one-frame baseband waveform (complex samples)
%   txData.txBits          — input bits before BCH+LDPC encoding
%   txData.txSymbols       — symbols at 1 sample/symbol (post matched filter)
%   txData.pilotValue      — known DVB-S2 pilot: (1+j)/sqrt(2)
%   txData.filterDelay     — matched filter group delay in symbols
%   txData.wavegen         — wavegen object (holds MODCOD/frame config)
% =========================================================================

%% Configure generator
wavegen = dvbs2WaveformGenerator;
wavegen.StreamFormat        = 'TS';
wavegen.NumInputStreams      = 1;
wavegen.FECFrame             = cfg.FECFrame;
wavegen.MODCOD               = cfg.CurrentMODCOD;
wavegen.HasPilots            = cfg.HasPilots;
wavegen.RolloffFactor        = cfg.RolloffFactor;
wavegen.FilterSpanInSymbols  = cfg.FilterSpanInSymbols;
wavegen.SamplesPerSymbol     = cfg.SamplesPerSymbol;

if isfield(cfg,'DFL') && ~isempty(cfg.DFL)
    wavegen.DFL = cfg.DFL;
end

% NumPLFrames = 1 — exactly one frame per call
NumPLFrames = 1;

%% Generate input bits (fresh random bits each call via PNSequence)
pn = comm.PNSequence( ...
    'Polynomial',         'x9+x5+1', ...
    'InitialConditions',  [zeros(1,8) 1], ...
    'VariableSizeOutput', true, ...
    'MaximumOutputSize',  [64800, 1]);

data    = cell(1, wavegen.NumInputStreams);
syncBits = [0;1;0;0;0;1;1;1];
pktLen   = 1496;

for i = 1:wavegen.NumInputStreams
    numPkts   = wavegen.MinNumPackets(i) * NumPLFrames;
    reset(pn);
    rawpkts   = pn(pktLen * numPkts);
    txRawPkts = reshape(rawpkts, pktLen, numPkts);
    txPkts    = [repmat(syncBits, 1, numPkts); txRawPkts];
    data{i}   = txPkts(:);
end

inputBits = data{1};

%% Generate waveform (BCH + LDPC + mapping + SRRC filtering)
in       = data;
waveform = [wavegen(in); flushFilter(wavegen)];

%% Downsample to symbol rate for pilot reference
txFilter = comm.RaisedCosineReceiveFilter( ...
    'RolloffFactor',         cfg.RolloffFactor, ...
    'FilterSpanInSymbols',   cfg.FilterSpanInSymbols, ...
    'InputSamplesPerSymbol', cfg.SamplesPerSymbol, ...
    'DecimationFactor',      cfg.SamplesPerSymbol);

txSym_raw   = txFilter(waveform);
filterDelay = cfg.FilterSpanInSymbols / 2;
txSymbols   = txSym_raw(filterDelay+1:end);
release(txFilter);

%% Pack output
txData.waveform    = waveform;
txData.txBits      = inputBits;
txData.txSymbols   = txSymbols;
txData.pilotValue  = (1/sqrt(2)) + 1j*(1/sqrt(2));
txData.filterDelay = filterDelay;
txData.wavegen     = wavegen;
