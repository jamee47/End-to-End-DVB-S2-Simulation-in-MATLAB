clc; clearvars; close all;

% Download LDPC parity matrices if not available
if ~exist('dvbs2xLDPCParityMatrices.mat','file')
    if ~exist('s2xLDPCParityMatrices.zip','file')
        url = 'https://ssd.mathworks.com/supportfiles/spc/satcom/DVB/s2xLDPCParityMatrices.zip';
        websave('s2xLDPCParityMatrices.zip', url);
        unzip('s2xLDPCParityMatrices.zip');
    end
    addpath('s2xLDPCParityMatrices');
end

% ── Es/No sweep ────────────────────────────────────────────────────────────
esno    = 8.8:0.01:9.05;       % <-- adjust range and step as needed
numEsNo = length(esno);

% ── Adaptive stopping thresholds ───────────────────────────────────────────
% Keep processing frames until minBitErrors is reached OR maxFrames is hit.
% This gives statistically reliable estimates at all BER levels:
%   high BER  → stops early (errors accumulate fast)
%   low BER   → runs longer (needs more frames to see errors)
minBitErrors = 200;     % minimum bit errors before accepting the estimate
maxFrames    = 1000;    % hard cap to prevent infinite loops at very low BER

berResults = nan(1, numEsNo);
perResults = nan(1, numEsNo);

% ── Transmitter config (set once — same for all Es/No) ────────────────────
wavegen = dvbs2WaveformGenerator;
wavegen.StreamFormat        = 'TS';
wavegen.NumInputStreams      = 1;
wavegen.FECFrame             = 'normal';
wavegen.MODCOD               = 18;          % 16APSK 2/3  ← change as needed
wavegen.HasPilots            = true;
wavegen.RolloffFactor        = 0.35;
wavegen.FilterSpanInSymbols  = 10;
wavegen.SamplesPerSymbol     = 2;
wavegen.DFL                  = getdfl(wavegen.MODCOD, wavegen.FECFrame);

% normFlag: per MathWorks reference — true only for APSK with "Outer radius as 1".
% dvbs2WaveformGenerator default ScalingMethod = "Unit average power" → false.
if wavegen.MODCOD >= 18
    normFlag = strcmpi(wavegen.ScalingMethod, 'Outer radius as 1');
else
    normFlag = false;
end

syncBits        = [0;1;0;0;0;1;1;1];
pktLen          = 1496;
pktTotal        = 8 + pktLen;                       % 1504 bits per TS packet
numPktsPerFrame = wavegen.MinNumPackets(1);
dataSize        = numPktsPerFrame * pktTotal;        % decoded bits per frame

% ── PHY frame geometry (fixed by MODCOD — compute once) ───────────────────
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

pilotLen    = numPilotBlks * pilotSymLen;
plFrameSize = dataLen + pilotLen + slotLen;

plScrambIntSeq = satcom.internal.dvbs.plScramblingIntegerSequence(0);
cMap  = [1; 1j; -1; -1j];
cSeq  = cMap(plScrambIntSeq + 1);

[~, pilotInd_raw] = satcom.internal.dvbs.pilotBlock(numPilotBlks);
pilotInd  = pilotInd_raw + slotLen;
refPilots = (1+1j)/sqrt(2) .* cSeq(pilotInd_raw);

filterDelay  = wavegen.FilterSpanInSymbols / 2;   % correct RRC delay (symbols)
searchRange  = 5;                                  % frame-sync search window

% ── Main Es/No loop ────────────────────────────────────────────────────────
for i = 1:numEsNo

    bitsErr       = 0;
    pktsErr       = 0;
    pktsRec       = 0;
    numFramesLost = 0;
    totalFrames   = 0;
    dataStInd     = 1;

    % Pre-generate a large bit pool — reused across the adaptive frame loop.
    % Pool size = maxFrames worth of packets.
    numPkts   = numPktsPerFrame * maxFrames;
    pn = comm.PNSequence('Polynomial','x9+x5+1', ...
        'InitialConditions', [zeros(1,8) 1], ...
        'VariableSizeOutput', true, ...
        'MaximumOutputSize', [pktLen * numPkts, 1]);
    reset(pn);
    rawpkts   = pn(pktLen * numPkts);
    txRawPkts = reshape(rawpkts, pktLen, numPkts);
    txPkts    = [repmat(syncBits, 1, numPkts); txRawPkts];
    inputBits = txPkts(:);

    snr = esno(i) - 10*log10(wavegen.SamplesPerSymbol);

    % ── Adaptive frame loop ────────────────────────────────────────────────
    % Process in batches of batchSize frames; regenerate waveform each batch.
    % Continue until minBitErrors is reached or maxFrames is exhausted.
    batchSize  = 20;
    batchStart = 1;    % which packet index the next batch starts from

    while bitsErr < minBitErrors && totalFrames < maxFrames

        % Number of frames to generate this batch
        framesThisBatch = min(batchSize, maxFrames - totalFrames);

        % Slice the pre-generated bit pool for this batch
        pktsThisBatch = numPktsPerFrame * framesThisBatch;
        pktOffset     = (batchStart - 1) * numPktsPerFrame;
        batchBits     = inputBits((pktOffset * pktTotal + 1) : ...
                                  (pktOffset + pktsThisBatch) * pktTotal);
        batchData     = {batchBits};

        % Generate waveform for this batch
        release(wavegen);
        wavegen.DFL = getdfl(wavegen.MODCOD, wavegen.FECFrame);
        txWaveform  = [wavegen(batchData); flushFilter(wavegen)];

        % AWGN channel
        rxWaveform = awgn(txWaveform, snr, 'measured');

        % Matched filter
        rxFilter = comm.RaisedCosineReceiveFilter( ...
            'RolloffFactor'        , wavegen.RolloffFactor, ...
            'FilterSpanInSymbols'  , wavegen.FilterSpanInSymbols, ...
            'InputSamplesPerSymbol', wavegen.SamplesPerSymbol, ...
            'DecimationFactor'     , wavegen.SamplesPerSymbol);
        rxSym_raw = rxFilter(rxWaveform);
        release(rxFilter);

        rxSymbols = rxSym_raw(filterDelay + 1 : end);
        rxSymbols = rxSymbols / sqrt(mean(abs(rxSymbols).^2));  % normalise

        % Process each frame in this batch
        for stIdx = 1 : framesThisBatch

            if bitsErr >= minBitErrors || totalFrames >= maxFrames
                break;
            end

            nominalStart = (stIdx - 1) * plFrameSize + 1;
            nominalEnd   = nominalStart + plFrameSize - 1;
            if nominalEnd + searchRange > length(rxSymbols)
                break;
            end

            % Coarse frame sync via pilot correlation
            bestOffset = 0;
            bestMetric = -Inf;
            for offset = -searchRange : searchRange
                cStart = nominalStart + offset;
                if cStart < 1 || (cStart + plFrameSize - 1) > length(rxSymbols)
                    continue;
                end
                candidate   = rxSymbols(cStart : cStart + plFrameSize - 1);
                pilotMetric = abs(mean(candidate(pilotInd) .* conj(refPilots)));
                if pilotMetric > bestMetric
                    bestMetric = pilotMetric;
                    bestOffset = offset;
                end
            end

            rxFrame = rxSymbols(nominalStart + bestOffset : ...
                                nominalStart + bestOffset + plFrameSize - 1);

            totalFrames = totalFrames + 1;

            % PL header recovery and validation
            phyParams = dvbsPLHeaderRecover(rxFrame(1:90), Mode="DVB-S2/S2X regular");
            M = phyParams.ModulationOrder;
            if M == 0; R = []; else; R = eval(phyParams.LDPCCodeIdentifier); end

            if M ~= modOrder || R ~= codeRate || ...
                    phyParams.FECFrameLength ~= cwLen || ~phyParams.HasPilots
                numFramesLost = numFramesLost + 1;
                continue;
            end

            % Noise variance estimate + decode
            nVar = DVBS2NoiseVarEstimate(rxFrame, pilotInd, refPilots, normFlag);
            [decBitsTemp, isFrameLost, pktCRC] = dvbs2BitRecover(rxFrame, nVar, normFlag);
            decBits = decBitsTemp{:};

            if ~isFrameLost && length(decBits) ~= dataSize
                isFrameLost = true;
            end

            numFramesLost = numFramesLost + isFrameLost;

            if ~isFrameLost
                % PER: CRC-based packet errors
                pktsErr = pktsErr + numel(pktCRC{:}) - sum(pktCRC{:});
                pktsRec = pktsRec + numel(pktCRC{:});

                % BER: per-frame bit comparison against correct portion of inputBits
                globalPktOffset = pktOffset + (stIdx - 1) * numPktsPerFrame;
                bitInd = globalPktOffset * pktTotal + (1 : dataSize);
                if bitInd(end) <= length(inputBits)
                    bitsErr = bitsErr + sum(inputBits(bitInd) ~= decBits);
                end
                dataStInd = dataStInd + 1;
            end
        end

        batchStart = batchStart + framesThisBatch;

        % If we've exhausted the pre-generated pool, stop
        if batchStart * numPktsPerFrame > numPkts
            break;
        end
    end

    % ── Aggregate BER and PER ──────────────────────────────────────────────
    if pktsRec > 0
        berResults(i) = bitsErr / (pktsRec * pktTotal);
        perResults(i) = pktsErr / pktsRec;
    else
        berResults(i) = 1.0;
        perResults(i) = 1.0;
    end

    % Flag low-confidence points (fewer errors than threshold)
    confFlag = '';
    if bitsErr < minBitErrors
        confFlag = '  [<200 errors — low confidence]';
    end

    fprintf('Es/No = %5.2f dB | Frames %d/%d lost | Errors %d bits | BER = %.3e | PER = %.3e%s\n', ...
        esno(i), numFramesLost, totalFrames, bitsErr, berResults(i), perResults(i), confFlag);

    clearvars -except esno numEsNo berResults perResults i ...
        wavegen modOrder codeRate cwLen dataLen slotLen pilotLen plFrameSize ...
        numPilotBlks pilotInd pilotInd_raw refPilots cSeq ...
        numPktsPerFrame dataSize pktTotal pktLen syncBits normFlag ...
        filterDelay searchRange minBitErrors maxFrames;
end

% ── Save ───────────────────────────────────────────────────────────────────
save('dvbs2_ber_results.mat', 'esno', 'berResults', 'perResults');
fprintf('\nSaved → dvbs2_ber_results.mat\n');

% ── Plot ───────────────────────────────────────────────────────────────────
validBER = berResults > 0;
validPER = perResults > 0;

figure;
semilogy(esno(validBER), berResults(validBER), 'b-o', ...
    'LineWidth', 1.5, 'MarkerFaceColor', 'b'); hold on;
semilogy(esno(validPER), perResults(validPER), 'r-s', ...
    'LineWidth', 1.5, 'MarkerFaceColor', 'r');
grid on;
legend('BER','PER','Location','southwest');
xlabel('E_s/N_0 (dB)');
ylabel('Error Rate');
title(sprintf('DVB-S2 BER & PER  (MODCOD %d, Normal FEC, Pilots ON)', wavegen.MODCOD));
ylim([1e-6 1]);

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