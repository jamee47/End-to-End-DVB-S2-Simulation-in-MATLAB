clc; clearvars; close all;

% ══════════════════════════════════════════════════════════════════════════
%  DVB-S2 Channel Estimation Dataset Generator
%  Saves to both:
%    dvbs2_channel_dataset.mat  — full precision struct array
%    dvbs2_channel_dataset.csv  — flat row-per-frame feature table
%
%  CSV column order (one row = one PL frame):
%    [plheader_re × 90 | data_re × Nd | pilot_re × Np |
%     plheader_im × 90 | data_im × Nd | pilot_im × Np |
%     H_LS_re × K | H_LS_im × K |
%     H_true_re | H_true_im |
%     snr_dB | nVar | modcod | Kfactor_dB]
%
%  NOTE: Nd (data symbols) ≈ 32 400 for QPSK normal frame.
%        CSV file size ≈ 0.5–1 GB per 2 000 frames at 6-digit precision.
%        Full raw symbols + per-symbol channel (H_true_sym) live in MAT only.
% ══════════════════════════════════════════════════════════════════════════

% ── Dataset control ────────────────────────────────────────────────────────
numFrames  = 10;        % ← total frames to generate
saveEvery  = 2;         % checkpoint save interval (MAT only)
matFile    = '../dataset_output/dvbs2_channel_dataset.mat';
csvFile    = '../dataset_output/dvbs2_channel_dataset.csv';

% ── SNR and channel diversity ──────────────────────────────────────────────
esno_pool      = [30];    % dB — drawn randomly per frame
specificRainAtt_dBpool = [0 3 5 7 10];         % Rician K (dB); 0≈Rayleigh, 10≈strong LOS

% ── Transmitter config ─────────────────────────────────────────────────────
wavegen = dvbs2WaveformGenerator;
wavegen.StreamFormat        = 'TS';
wavegen.NumInputStreams      = 1;
wavegen.FECFrame             = 'normal';
wavegen.MODCOD               = 5;       % ← QPSK 3/5; change as needed
wavegen.HasPilots            = true;
wavegen.RolloffFactor        = 0.35;
wavegen.FilterSpanInSymbols  = 10;
wavegen.SamplesPerSymbol     = 2;
wavegen.DFL                  = getdfl(wavegen.MODCOD, wavegen.FECFrame);

if wavegen.MODCOD >= 18
    normFlag = strcmpi(wavegen.ScalingMethod, 'Outer radius as 1');
else
    normFlag = false;
end

syncBits        = [0;1;0;0;0;1;1;1];
pktLen          = 1496;
pktTotal        = 8 + pktLen;
numPktsPerFrame = wavegen.MinNumPackets(1);

% ── PHY frame geometry ─────────────────────────────────────────────────────
[modOrder, codeRate, cwLen] = satcom.internal.dvbs.getS2PHYParams( ...
    wavegen.MODCOD, wavegen.FECFrame);

dataLen      = cwLen / log2(modOrder);
slotLen      = 90;
pilotBlkFreq = 16;
pilotSymLen  = 36;

numPilotBlks = floor(dataLen / (slotLen * pilotBlkFreq));
if floor(dataLen/(slotLen*16)) == dataLen/(slotLen*pilotBlkFreq)
    numPilotBlks = numPilotBlks - 1;
end

pilotLen    = numPilotBlks * pilotSymLen;    % total pilot symbols
plFrameSize = dataLen + pilotLen + slotLen;  % total symbols per PL frame

% PL scrambling
plScrambIntSeq = satcom.internal.dvbs.plScramblingIntegerSequence(0);
cMap  = [1; 1j; -1; -1j];
cSeq  = cMap(plScrambIntSeq + 1);

% Pilot indices and references — force to row vectors to avoid dim mismatch
[~, pilotInd_raw] = satcom.internal.dvbs.pilotBlock(numPilotBlks);
pilotInd_raw = pilotInd_raw(:)';                              % [1 × total_pilot_syms]
pilotInd     = pilotInd_raw + slotLen;                        % row, within full frame
refPilots    = (1+1j)/sqrt(2) .* cSeq(pilotInd_raw(:));      % col [total_pilot_syms × 1]
refPilotMat  = reshape(refPilots, pilotSymLen, numPilotBlks).'; % [numPilotBlks × 36]

% Index maps
headerInd = 1 : slotLen;                                              % [1 × 90]
dataInd   = setdiff(1:plFrameSize, [headerInd(:)', pilotInd(:)']);    % [1 × Nd]

Nd = length(dataInd);           % number of data symbols per frame
Np = length(pilotInd);          % total pilot symbols per frame  (K × 36)
K  = numPilotBlks;              % number of pilot blocks

filterDelay = wavegen.FilterSpanInSymbols / 2;
sampleRate  = wavegen.SamplesPerSymbol;

% ── Metadata struct ────────────────────────────────────────────────────────
metadata.plFrameSize     = plFrameSize;
metadata.numPilotBlks    = numPilotBlks;
metadata.pilotSymLen     = pilotSymLen;
metadata.pilotInd        = pilotInd;
metadata.pilotInd_raw    = pilotInd_raw;
metadata.refPilots       = refPilots;
metadata.refPilotMat     = refPilotMat;
metadata.dataInd         = dataInd;
metadata.headerInd       = headerInd;
metadata.slotLen         = slotLen;
metadata.modcod          = wavegen.MODCOD;
metadata.modOrder        = modOrder;
metadata.codeRate        = codeRate;
metadata.cwLen           = cwLen;
metadata.dataLen         = dataLen;
metadata.Nd              = Nd;
metadata.Np              = Np;
metadata.numPktsPerFrame = numPktsPerFrame;
metadata.pktTotal        = pktTotal;

% ── Pre-allocate MAT dataset struct ───────────────────────────────────────
dataset(numFrames) = struct( ...
    'rxFrame',    [], ...   % [plFrameSize × 1] complex — full frame
    'rxDataSyms', [], ...   % [Nd × 1] complex — data symbols
    'rxPilotMat', [], ...   % [K × 36] complex — pilot matrix
    'H_LS',       [], ...   % [K × 1] complex — LS estimate per block
    'H_true',     [], ...   % complex scalar — mean channel (BLSTM label)
    'H_true_sym', [], ...   % [plFrameSize × 1] complex — per-symbol (Transformer label)
    'nVar',       [], ...
    'esno_dB',    [], ...
    'modcod',     [], ...
    'rainAtt_dB', []);

% ── Open CSV and write header row ─────────────────────────────────────────
% Column names mirror the stored data exactly.
fid = fopen(csvFile, 'w');
if fid == -1
    error('Cannot open CSV file for writing: %s', csvFile);
end

% Build header string
hdr = '';
for c = 1:90;  hdr = [hdr, sprintf('plheader_re_%d,', c)]; end
for c = 1:Nd;  hdr = [hdr, sprintf('data_re_%d,', c)];    end
for c = 1:Np;  hdr = [hdr, sprintf('pilot_re_%d,', c)];   end
for c = 1:90;  hdr = [hdr, sprintf('plheader_im_%d,', c)]; end
for c = 1:Nd;  hdr = [hdr, sprintf('data_im_%d,', c)];    end
for c = 1:Np;  hdr = [hdr, sprintf('pilot_im_%d,', c)];   end
for c = 1:K;   hdr = [hdr, sprintf('H_LS_re_%d,', c)];    end
for c = 1:K;   hdr = [hdr, sprintf('H_LS_im_%d,', c)];    end
hdr = [hdr, 'H_true_re,H_true_im,snr_dB,nVar,modcod,rainAtt_dB'];
fprintf(fid, '%s\n', hdr);

fprintf('CSV header written — %d columns per row.\n', ...
    90 + Nd + Np + 90 + Nd + Np + K + K + 6);
fprintf('Starting dataset generation: %d frames, MODCOD %d\n\n', ...
    numFrames, wavegen.MODCOD);
fprintf('%-6s %-10s %-10s %-12s %-12s\n', ...
    'Frame','Es/No(dB)','rainAtt(dB)','|H_true|','nVar');
fprintf('%s\n', repmat('-',1,52));

% ── Main generation loop ───────────────────────────────────────────────────
for frameIdx = 1 : numFrames

    % Random channel conditions for this frame
    esno_dB = esno_pool(randi(length(esno_pool)));
    rainAtt_dB = specificRainAtt_dBpool(randi(length(specificRainAtt_dBpool)));
    snr     = esno_dB - 10*log10(sampleRate);

    % ── Input bits ────────────────────────────────────────────────────────
    pn = comm.PNSequence('Polynomial','x9+x5+1', ...
        'InitialConditions', [zeros(1,8) 1], ...
        'VariableSizeOutput', true, ...
        'MaximumOutputSize', [pktLen * numPktsPerFrame, 1]);
    reset(pn);
    rawBits   = pn(pktLen * numPktsPerFrame);
    txRawPkts = reshape(rawBits, pktLen, numPktsPerFrame);
    txPkts    = [repmat(syncBits, 1, numPktsPerFrame); txRawPkts];
    inputBits = txPkts(:);

    % ── TX waveform ───────────────────────────────────────────────────────
    release(wavegen);
    wavegen.DFL = getdfl(wavegen.MODCOD, wavegen.FECFrame);
    txWaveform  = [wavegen({inputBits}); flushFilter(wavegen)];

    % % ── Rician fading channel ──────────────────────────────────────────────
    % Ap_linear = 10^(-Ap_dB/20);
    % phi = 2*pi*rand(1);
    %
    % h_tilde = sqrt(Ap_linear)*exp(-1i*phi);
    % attWaveform = h_tilde*txWaveform;
    %
    % snr = esno-10*log10(SPS);
    %
    % fadedWaveform = awgn(attWaveform,snr,'measured');

    rainAtt = 10^(-rainAtt_dB/20);
    phi = 2*pi*rand(1);

    h_tilde = sqrt(rainAtt)*exp(1i*phi);
    rxFaded = h_tilde*txWaveform;

    % ── AWGN ──────────────────────────────────────────────────────────────
    rxNoisy = awgn(rxFaded, snr, 'measured');
    clearvars txWaveform rxFaded rawBits txRawPkts txPkts inputBits pn ricianChan;

    % ── Matched filter ─────────────────────────────────────────────────────
    rxFilter = comm.RaisedCosineReceiveFilter( ...
        'RolloffFactor'        , wavegen.RolloffFactor, ...
        'FilterSpanInSymbols'  , wavegen.FilterSpanInSymbols, ...
        'InputSamplesPerSymbol', sampleRate, ...
        'DecimationFactor'     , sampleRate);
    rxSym_raw = rxFilter(rxNoisy);
    release(rxFilter);
    clearvars rxNoisy rxFilter;

    % Discard transient and normalise
    rxSymbols = rxSym_raw(filterDelay + 1 : end);
    clearvars rxSym_raw;
    rxSymbols = rxSymbols / sqrt(mean(abs(rxSymbols).^2));

    % Per-symbol channel aligned to rxSymbols
    %gainsSym  = pathGains(sampleRate : sampleRate : end);
    %gainsSym  = gainsSym(filterDelay + 1 : end);
    %clearvars pathGains;

    % if length(rxSymbols) < plFrameSize || length(gainsSym) < plFrameSize
    %     warning('Frame %d: insufficient symbols — skipping.', frameIdx);
    %     clearvars rxSymbols gainsSym;
    %     continue;
    % end

    % ── Extract frame components ───────────────────────────────────────────
    rxFrame     = rxSymbols(1 : plFrameSize);               % full PL frame
    %H_true_sym  = gainsSym(1 : plFrameSize);                % per-symbol channel
    %clearvars rxSymbols gainsSym;

    rxHeader    = rxFrame(headerInd);                        % [90 × 1]
    rxDataSyms  = rxFrame(dataInd);                          % [Nd × 1]
    rxPilotVec  = rxFrame(pilotInd);                         % [Np × 1] (row→col ok)
    rxPilotMat  = reshape(rxFrame(pilotInd(:)), ...
        pilotSymLen, numPilotBlks).';      % [K × 36]

    % LS estimate per pilot block (averaged over 36 symbols)
    H_LS = mean(rxPilotMat .* conj(refPilotMat), 2);        % [K × 1]

    % True channel scalar
    H_true = h_tilde;

    % Noise variance
    nVar = DVBS2NoiseVarEstimate(rxFrame, pilotInd, refPilots, normFlag);

    % ── Store in MAT dataset ──────────────────────────────────────────────
    dataset(frameIdx).rxFrame      = rxFrame;
    dataset(frameIdx).rxDataSyms   = rxDataSyms;
    dataset(frameIdx).rxPilotMat   = rxPilotMat;
    dataset(frameIdx).H_LS         = H_LS;
    dataset(frameIdx).H_true       = h_tilde;
    dataset(frameIdx).H_true_sym   = h_tilde;
    dataset(frameIdx).nVar         = nVar;
    dataset(frameIdx).esno_dB      = esno_dB;
    dataset(frameIdx).modcod       = wavegen.MODCOD;
    dataset(frameIdx).Kfactor_dB   = 10*log10(rainAtt_dB);

    % ── Write one CSV row ─────────────────────────────────────────────────
    % Column order:
    %   plheader_re | data_re | pilot_re |
    %   plheader_im | data_im | pilot_im |
    %   H_LS_re | H_LS_im |
    %   H_true_re | H_true_im | snr_dB | nVar | modcod | Kfactor_dB
    rowVec = [ ...
        real(rxHeader(:))',  real(rxDataSyms(:))',  real(rxPilotVec(:))', ...
        imag(rxHeader(:))',  imag(rxDataSyms(:))',  imag(rxPilotVec(:))', ...
        real(H_LS(:))',      imag(H_LS(:))', ...
        real(H_true),        imag(H_true), ...
        esno_dB,             nVar, ...
        wavegen.MODCOD,      rainAtt_dB ...
        ];

    % Write row: values separated by commas, 6 decimal places
    fprintf(fid, '%s\n', strjoin(arrayfun(@(v) sprintf('%.6g', v), ...
        rowVec, 'UniformOutput', false), ','));

    fprintf('%-6d %-10.1f %-10.1f %-12.4f %-12.4e\n', ...
        frameIdx, esno_dB, 10*log10(rainAtt_dB), abs(H_true), nVar);

    % ── Clear per-frame temporaries ───────────────────────────────────────
    clearvars rxFrame rxHeader rxDataSyms rxPilotVec rxPilotMat ...
        H_LS H_true H_true_sym nVar rowVec esno_dB Kfactor snr;

    % ── MAT checkpoint ────────────────────────────────────────────────────
    if mod(frameIdx, saveEvery) == 0 || frameIdx == numFrames
        metadata.framesGenerated = frameIdx;
        save(matFile, 'dataset', 'metadata', '-v7.3');
        fprintf('  → MAT checkpoint: %d/%d frames → %s\n', ...
            frameIdx, numFrames, matFile);
    end
end

fclose(fid);

fprintf('\nDataset generation complete.\n');
fprintf('  MAT : %s\n', matFile);
fprintf('  CSV : %s  (%d cols × %d rows)\n', csvFile, ...
    90 + Nd + Np + 90 + Nd + Np + K + K + 6, numFrames);

% ══════════════════════════════════════════════════════════════════════════
%  Helper functions
% ══════════════════════════════════════════════════════════════════════════

function dfl = getdfl(modCod, fecFrame)
if strcmp(fecFrame, 'normal')
    nDefVal = [16008 21408 25728 32208 38688 43040 48408 51648 53840 57472 ...
        58192 38688 43040 48408 53840 57472 58192 43040 48408 51648 ...
        53840 57472 58192 48408 51648 53840 57472 58192] - 80;
else
    nDefVal = [3072 5232 6312 7032 9552 10632 11712 12432 13152 14232 0 ...
        9552 10632 11712 13152 14232 0 10632 11712 12432 13152 14232 ...
        0 11712 12432 13152 14232 0] - 80;
end
dfl = nDefVal(modCod);
end

function nVarEst = DVBS2NoiseVarEstimate(rxData, pilotInd, refPilots, normFlag)
if normFlag
    rxData = rxData / sqrt(mean(abs(rxData).^2));
end
rxPilots = rxData(pilotInd(:));
pSigN    = mean(abs(rxPilots).^2);
pSig     = abs(mean(rxPilots .* conj(refPilots(:)))).^2;
nVarEst  = abs(pSigN - pSig);
end