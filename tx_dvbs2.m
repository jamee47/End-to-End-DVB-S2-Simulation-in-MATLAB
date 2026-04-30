function txData = tx_dvbs2(cfg)
% =========================================================================
% TX_DVBS2 — DVB-S2 Transmitter Module
% =========================================================================
% Generates a DVB-S2 waveform for the given configuration.
% Returns a struct containing the waveform, transmitted symbols,
% pilot positions, and original bits for BER computation.
%
% INPUT:
%   cfg  — configuration struct from main_generate_dataset.m
%
% OUTPUT:
%   txData.waveform        — filtered, modulated waveform (complex samples)
%   txData.txBits          — original information bits (for BER)
%   txData.txSymbols       — complex symbols after downsampling (reference)
%   txData.pilotValue      — known pilot symbol value (1/sqrt(2) + j/sqrt(2))
%   txData.wavegen         — dvbs2WaveformGenerator object (kept for reference)
%   txData.NumFrames       — number of PL frames generated
%   txData.SymbolsPerFrame — estimated symbols per PLFRAME
% =========================================================================

%% Configure DVB-S2 Waveform Generator
wavegen = dvbs2WaveformGenerator;
wavegen.StreamFormat         = 'TS';
wavegen.NumInputStreams       = 1;
wavegen.FECFrame              = cfg.FECFrame;
wavegen.MODCOD                = cfg.CurrentMODCOD;
wavegen.HasPilots             = cfg.HasPilots;
wavegen.RolloffFactor         = cfg.RolloffFactor;
wavegen.FilterSpanInSymbols   = cfg.FilterSpanInSymbols;
wavegen.SamplesPerSymbol      = cfg.SamplesPerSymbol;

% DFL: Data Field Length — use maximum for short frame
% Short frame FECFRAME = 16200 bits
% Max DFL for short frame TS mode
wavegen.DFL = 3008;   % valid DFL for short frame TS

%% Generate Input Bits (Transport Stream packets)
syncBits = [0;1;0;0;0;1;1;1];
pktLen   = 1496;

pn = comm.PNSequence('Polynomial','x9+x5+1', ...
    'InitialConditions',[zeros(1,8) 1], ...
    'VariableSizeOutput', true, ...
    'MaximumOutputSize', [64800 * cfg.NumFrames, 1]);

numPkts   = wavegen.MinNumPackets(1) * cfg.NumFrames;
reset(pn);
rawpkts   = pn(pktLen * numPkts);
txRawPkts = reshape(rawpkts, pktLen, numPkts);
txPkts    = [repmat(syncBits, 1, numPkts); txRawPkts];
inputBits = txPkts(:);

data = {inputBits};

%% Generate Waveform
waveform = [wavegen(data); flushFilter(wavegen)];

%% Extract TX Symbols (downsample for pilot reference)
% Apply matched filter to get symbol-rate samples
txFilter = comm.RaisedCosineReceiveFilter( ...
    'RolloffFactor',         cfg.RolloffFactor, ...
    'FilterSpanInSymbols',   cfg.FilterSpanInSymbols, ...
    'InputSamplesPerSymbol', cfg.SamplesPerSymbol, ...
    'DecimationFactor',      cfg.SamplesPerSymbol);

txSym_raw    = txFilter(waveform);
filterDelay  = cfg.FilterSpanInSymbols / 2;
txSymbols    = txSym_raw(filterDelay+1:end);
release(txFilter);

%% Known Pilot Symbol Value (DVB-S2 standard)
% Pilots are un-modulated: I = 1/sqrt(2), Q = 1/sqrt(2)
pilotValue = (1/sqrt(2)) + 1j*(1/sqrt(2));

%% Estimate symbols per PLFRAME for reference
% Short frame: 16200 / bitsPerSymbol + overhead
totalSymbols    = length(txSymbols);
symbolsPerFrame = floor(totalSymbols / cfg.NumFrames);

%% Pack output
txData.waveform        = waveform;
txData.txBits          = inputBits;
txData.txSymbols       = txSymbols;
txData.pilotValue      = pilotValue;
txData.wavegen         = wavegen;
txData.NumFrames       = cfg.NumFrames;
txData.SymbolsPerFrame = symbolsPerFrame;
txData.filterDelay     = filterDelay;

fprintf('    [TX] Waveform length : %d samples\n', length(waveform));
fprintf('    [TX] Total TX symbols: %d\n', totalSymbols);
fprintf('    [TX] ~Symbols/frame  : %d\n', symbolsPerFrame);
% Fs = 1e+06; 
% spectrum = spectrumAnalyzer(SampleRate=Fs, AveragingMethod='exponential', ForgettingFactor=1);
% spectrum(waveform);
% release(spectrum);
% 
% % Constellation Diagram
% constel = comm.ConstellationDiagram('ColorFading', true, ...
%     'ShowTrajectory', 0, ...
%     'SamplesPerSymbol', 4, ...
%     'ShowReferenceConstellation', false);    
% 
% constel(waveform);
% release(constel);
end