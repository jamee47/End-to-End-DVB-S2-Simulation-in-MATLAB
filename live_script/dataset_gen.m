clc; clearvars; close all;

% ── Dataset control ────────────────────────────────────────────────────────
numFrames  = 1;
saveEvery  = 1;
outputFile = 'dvbs2_channel_dataset.csv';

% ── SNR and channel diversity ──────────────────────────────────────────────
esno_pool      = [5 10 15 20 25 30];
Kfactor_dBpool = [0 3 5 7 10];

% ── Transmitter config ─────────────────────────────────────────────────────
wavegen = dvbs2WaveformGenerator;
wavegen.StreamFormat        = 'TS';
wavegen.NumInputStreams      = 1;
wavegen.FECFrame             = 'normal';
wavegen.MODCOD               = 5;
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
numPktsPerFrame = wavegen.MinNumPackets(1);

[modOrder, codeRate, cwLen] = satcom.internal.dvbs.getS2PHYParams( ...
    wavegen.MODCOD, wavegen.FECFrame);

dataLen      = cwLen / log2(modOrder);
slotLen      = 90;
pilotSymLen  = 36;

numPilotBlks = floor(dataLen / (slotLen * 16));
if floor(dataLen/(slotLen*16)) == dataLen/(slotLen*16)
    numPilotBlks = numPilotBlks - 1;
end

pilotLen    = numPilotBlks * pilotSymLen;
plFrameSize = dataLen + pilotLen + slotLen;

plScrambIntSeq = satcom.internal.dvbs.plScramblingIntegerSequence(0);
cMap  = [1; 1j; -1; -1j];
cSeq  = cMap(plScrambIntSeq + 1);

[~, pilotInd_raw] = satcom.internal.dvbs.pilotBlock(numPilotBlks);
pilotInd_raw = pilotInd_raw(:)';
pilotInd     = pilotInd_raw + slotLen;
refPilots    = (1+1j)/sqrt(2) .* cSeq(pilotInd_raw(:));
refPilotMat  = reshape(refPilots, pilotSymLen, numPilotBlks).';

headerInd = 1 : slotLen;
dataInd   = setdiff(1:plFrameSize, [headerInd(:)', pilotInd(:)']);
nData     = length(dataInd);

filterDelay = wavegen.FilterSpanInSymbols / 2;
sampleRate  = wavegen.SamplesPerSymbol;

% ════════════════════════════════════════════════════════════════════════
%  Build header then open file
%  Naming convention:
%    scalars          → plain column name  e.g. esno_dB
%    complex scalar   → H_true_re, H_true_im
%    complex vector   → H_LS_re_0 .. H_LS_re_N, H_LS_im_0 .. H_LS_im_N
%    complex matrix   → rxPilotMat_re_b0_s0 ... (block b, symbol s)
% ════════════════════════════════════════════════════════════════════════

fid = fopen(outputFile, 'w');
if fid == -1, error('Cannot open output file: %s', outputFile); end

% --- scalar columns ---
header = 'frame_idx,esno_dB,Kfactor_dB,modcod,nVar,H_true_re,H_true_im';

% --- H_LS: numPilotBlks complex values ---
for k = 0:numPilotBlks-1
    header = [header, sprintf(',H_LS_re_%d', k)];
end
for k = 0:numPilotBlks-1
    header = [header, sprintf(',H_LS_im_%d', k)];
end

% --- H_true_sym: plFrameSize complex values ---
for k = 0:plFrameSize-1
    header = [header, sprintf(',H_true_sym_re_%d', k)];
end
for k = 0:plFrameSize-1
    header = [header, sprintf(',H_true_sym_im_%d', k)];
end

% --- rxFrame: plFrameSize complex values ---
for k = 0:plFrameSize-1
    header = [header, sprintf(',rxFrame_re_%d', k)];
end
for k = 0:plFrameSize-1
    header = [header, sprintf(',rxFrame_im_%d', k)];
end

% --- rxDataSyms: nData complex values ---
for k = 0:nData-1
    header = [header, sprintf(',rxDataSyms_re_%d', k)];
end
for k = 0:nData-1
    header = [header, sprintf(',rxDataSyms_im_%d', k)];
end

% --- rxPilotMat: numPilotBlks x pilotSymLen complex values ---
for b = 0:numPilotBlks-1
    for s = 0:pilotSymLen-1
        header = [header, sprintf(',rxPilotMat_re_b%d_s%d', b, s)];
    end
end
for b = 0:numPilotBlks-1
    for s = 0:pilotSymLen-1
        header = [header, sprintf(',rxPilotMat_im_b%d_s%d', b, s)];
    end
end

fprintf(fid, '%s\n', header);

fprintf('Starting dataset generation: %d frames, MODCOD %d\n', numFrames, wavegen.MODCOD);
fprintf('%-6s %-10s %-10s %-12s %-12s\n','Frame','Es/No(dB)','K(dB)','|H_true|','nVar');
fprintf('%s\n', repmat('-',1,52));

% ── Main generation loop ──────────────────────────────────────────────────
for frameIdx = 1 : numFrames

    esno_dB = esno_pool(randi(length(esno_pool)));
    Kfactor = 10^(Kfactor_dBpool(randi(length(Kfactor_dBpool)))/10);
    snr     = esno_dB - 10*log10(sampleRate);

    pn = comm.PNSequence('Polynomial','x9+x5+1', ...
        'InitialConditions',[zeros(1,8) 1], ...
        'VariableSizeOutput',true, ...
        'MaximumOutputSize',[pktLen * numPktsPerFrame, 1]);
    reset(pn);
    rawBits   = pn(pktLen * numPktsPerFrame);
    txRawPkts = reshape(rawBits, pktLen, numPktsPerFrame);
    txPkts    = [repmat(syncBits,1,numPktsPerFrame); txRawPkts];
    inputBits = txPkts(:);

    release(wavegen);
    wavegen.DFL = getdfl(wavegen.MODCOD, wavegen.FECFrame);
    txWaveform  = [wavegen({inputBits}); flushFilter(wavegen)];

    ricianChan = comm.RicianChannel( ...
        'SampleRate',          sampleRate, ...
        'PathDelays',          0, ...
        'AveragePathGains',    0, ...
        'KFactor',             Kfactor, ...
        'MaximumDopplerShift', 0.01, ...
        'PathGainsOutputPort', true);

    [rxFaded, pathGains] = ricianChan(txWaveform);
    rxNoisy = awgn(rxFaded, snr, 'measured');
    clearvars txWaveform rxFaded pn rawBits txRawPkts txPkts inputBits ricianChan;

    rxFilter = comm.RaisedCosineReceiveFilter( ...
        'RolloffFactor',         wavegen.RolloffFactor, ...
        'FilterSpanInSymbols',   wavegen.FilterSpanInSymbols, ...
        'InputSamplesPerSymbol', sampleRate, ...
        'DecimationFactor',      sampleRate);
    rxSym_raw = rxFilter(rxNoisy);
    release(rxFilter);
    clearvars rxNoisy rxFilter;

    rxSymbols = rxSym_raw(filterDelay + 1 : end);
    clearvars rxSym_raw;
    rxSymbols = rxSymbols / sqrt(mean(abs(rxSymbols).^2));

    gainsSym = pathGains(sampleRate : sampleRate : end);
    gainsSym = gainsSym(filterDelay + 1 : end);
    clearvars pathGains;

    if length(rxSymbols) < plFrameSize || length(gainsSym) < plFrameSize
        warning('Frame %d: insufficient symbols — skipping.', frameIdx);
        clearvars rxSymbols gainsSym;
        continue;
    end

    rxFrame    = rxSymbols(1 : plFrameSize);
    H_true_sym = gainsSym(1 : plFrameSize);
    clearvars rxSymbols gainsSym;

    rxDataSyms = rxFrame(dataInd);
    rxPilotMat = reshape(rxFrame(pilotInd), pilotSymLen, numPilotBlks).';
    H_LS       = mean(rxPilotMat .* conj(refPilotMat), 2);
    H_true     = mean(H_true_sym);
    nVar       = DVBS2NoiseVarEstimate(rxFrame, pilotInd, refPilots, normFlag);
    Kfactor_dB = 10*log10(Kfactor);

    % ── Assemble one row: all reals first within each field, then imags ───
    row = [ frameIdx, esno_dB, Kfactor_dB, double(wavegen.MODCOD), nVar, ...
            real(H_true), imag(H_true), ...
            real(H_LS(:))',    imag(H_LS(:))', ...
            real(H_true_sym(:))', imag(H_true_sym(:))', ...
            real(rxFrame(:))',    imag(rxFrame(:))', ...
            real(rxDataSyms(:))', imag(rxDataSyms(:))', ...
            real(rxPilotMat(:))', imag(rxPilotMat(:))' ];

    % Write row — use %.8g for compact but precise float representation
    fprintf(fid, '%.8g', row(1));
    fprintf(fid, ',%.8g', row(2:end));
    fprintf(fid, '\n');

    fprintf('%-6d %-10.1f %-10.1f %-12.4f %-12.4e\n', ...
        frameIdx, esno_dB, Kfactor_dB, abs(H_true), nVar);

    clearvars rxFrame rxDataSyms rxPilotMat H_LS H_true H_true_sym ...
              nVar esno_dB Kfactor Kfactor_dB snr row;

    if mod(frameIdx, saveEvery) == 0 || frameIdx == numFrames
        fprintf('  → Checkpoint: %d/%d frames written to %s\n', ...
            frameIdx, numFrames, outputFile);
    end
end

fclose(fid);
fprintf('\nDataset generation complete. Output: %s\n', outputFile);

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
    rxPilots = rxData(pilotInd, :);
    pSigN    = mean(abs(rxPilots).^2);
    pSig     = abs(mean(rxPilots .* conj(refPilots))).^2;
    nVarEst  = abs(pSigN - pSig);
end