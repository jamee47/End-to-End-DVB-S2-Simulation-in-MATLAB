%% ber_analysis_ls.m
%  DVB-S2 BER/PER sweep over ITU-R P.618 tropical rain fading channel
%  with LS pilot-aided channel estimation and equalization.
%
%  Key differences from ber_analysis.m:
%   1. Channel: rain fading (h_k, Eq.7–11 of paper) + AWGN, replacing
%      plain AWGN-only channel.
%   2. AWGN noise is calibrated to the post-fade received power so that
%      Es/N0 refers to the signal level at the receiver, not the
%      transmitter (critical under heavy tropical attenuation).
%   3. LS channel estimation (Eq.12–16) is applied per frame before
%      calling dvbs2BitRecover.  The raw frame is NOT pre-normalised;
%      equalization is done by dividing by H_ls.
%   4. NMSE is tracked alongside BER/PER.
%   5. The pilot correlation frame-sync metric is run on the
%      un-equalized frame; equalization happens afterwards, just before
%      decoding — matching the paper's receiver chain order.
%   6. A separate Es/No axis is printed showing the *effective*
%      post-fade Es/N0 so you can compare with paper Fig. 8.
%
%  Channel scenario: tropical (Penang, Malaysia) at p = 0.01%.
%  To switch to temperate (Athens), change chanParams below.

clc; clearvars; close all;

%% =========================================================
%  LDPC PARITY MATRICES
%% =========================================================
if ~exist('dvbs2xLDPCParityMatrices.mat','file')
    if ~exist('s2xLDPCParityMatrices.zip','file')
        url = 'https://ssd.mathworks.com/supportfiles/spc/satcom/DVB/s2xLDPCParityMatrices.zip';
        websave('s2xLDPCParityMatrices.zip', url);
        unzip('s2xLDPCParityMatrices.zip');
    end
    addpath('s2xLDPCParityMatrices');
end

%% =========================================================
%  CHANNEL PARAMETERS — Tropical (Penang, Malaysia)
%  To simulate temperate channel, replace with Athens values:
%    chanParams.Latitude          = 37.90
%    chanParams.Altitude_m        = 15
%    chanParams.R001_mmph         = 24
%    chanParams.IsothermHeight_km = 3.0
%% =========================================================
chanParams.Frequency_GHz     = 20;    % Ka-band, GHz
chanParams.Latitude          = 5.17;  % deg N
chanParams.Altitude_m        = 57;    % m
chanParams.ElevationAngle    = 45;    % deg
chanParams.Polarization      = 'V';
chanParams.R001_mmph         = 130;   % mm/h  (tropical)
chanParams.IsothermHeight_km = 4.5;   % km    (tropical)
pRain                        = 0.01;  % time-percentage for exceedance

%% =========================================================
%  SATELLITE LINK BUDGET PARAMETERS
%% =========================================================
sysP.G_s_dBi        = 52.0;
sysP.G_r_dBi        = 41.7;
sysP.NoiseBW_Hz     = 50e6;
sysP.NoiseTemp_K    = 207;
sysP.Theta3dB_deg   = 0.4;
sysP.BeamOffset_deg = 0.0;

%% =========================================================
%  COMPUTE RAIN FADING CHANNEL ONCE
%  The channel is fixed for the entire sweep (slow-fading / per-run
%  realisation).  Re-run the script for Monte-Carlo averaging.
%% =========================================================
Ap_dB   = calcRainAttenuationP618(chanParams, pRain);
Ap      = 10^(Ap_dB / 20);          % amplitude attenuation factor
phi     = 2*pi*rand(1);              % random carrier phase

% Rain fading coefficient (Eq. 8): attenuation direction → 1/Ap
h_rain  = (1/Ap^0.5) * exp(-1j*phi);

% Free-space + receive-antenna factor (Eq. 9)
c_light = 3e8;
kb      = 1.38064852e-23;
freq_hz = chanParams.Frequency_GHz * 1e9;
lambda  = c_light / freq_hz;
d0      = 35788e3;
Gr      = 10^(sysP.G_r_dBi / 10);
b_max   = (lambda/(4*pi))^2 * (1/(d0^2)) * ...
          (Gr / (kb * sysP.NoiseBW_Hz * sysP.NoiseTemp_K));

% Satellite beam gain (Eq. 10) — boresight → pattern = 1
if sysP.BeamOffset_deg < 1e-9
    b_gain = 10^(sysP.G_s_dBi / 10);
else
    Gs      = 10^(sysP.G_s_dBi / 10);
    theta   = deg2rad(sysP.BeamOffset_deg);
    t3dB    = deg2rad(sysP.Theta3dB_deg);
    u_val   = 2.07123 * sin(theta) / sin(t3dB);
    bPat    = besselj(1,u_val)/(2*u_val) + 36*besselj(3,u_val)/(u_val^3);
    b_gain  = Gs * abs(bPat)^2;
end

% Full channel coefficient (Eq. 11)
h_k     = h_rain * sqrt(b_gain) * sqrt(b_max);
chanPow_dB = 20*log10(abs(h_k));

fprintf('========== CHANNEL SUMMARY ==========\n');
fprintf('Rain attenuation Ap_dB : %.4f dB\n', Ap_dB);
fprintf('Channel coeff |h_k|    : %.4e\n',    abs(h_k));
fprintf('Channel power          : %.2f dB\n', chanPow_dB);
fprintf('=====================================\n\n');

%% =========================================================
%  Es/No SWEEP PARAMETERS
%  Note: esno here is the *requested* Es/N0 at the receiver
%  input (post-fade).  The noise power is set relative to
%  the received signal power, so this is what the decoder sees.
%% =========================================================
esno         = 8.8 : 0.01 : 9.05;   % dB  — adjust range/step as needed
numEsNo      = length(esno);
minBitErrors = 200;
maxFrames    = 1000;

berResults  = nan(1, numEsNo);
perResults  = nan(1, numEsNo);
nmseResults = nan(1, numEsNo);

%% =========================================================
%  TRANSMITTER CONFIG (fixed for all Es/No points)
%% =========================================================
wavegen = dvbs2WaveformGenerator;
wavegen.StreamFormat       = 'TS';
wavegen.NumInputStreams     = 1;
wavegen.FECFrame            = 'normal';
wavegen.MODCOD              = 18;      % 16APSK 2/3 — change as needed
wavegen.HasPilots           = true;
wavegen.RolloffFactor       = 0.35;
wavegen.FilterSpanInSymbols = 10;
wavegen.SamplesPerSymbol    = 2;
wavegen.DFL                 = getdfl(wavegen.MODCOD, wavegen.FECFrame);

% normFlag controls power normalisation passed to dvbs2BitRecover.
% For LS equalization the frame is already equalised; we still pass
% normFlag to match the decoder's expected scaling convention.
if wavegen.MODCOD >= 18
    normFlag = strcmpi(wavegen.ScalingMethod, 'Outer radius as 1');
else
    normFlag = false;
end

syncBits        = [0;1;0;0;0;1;1;1];
pktLen          = 1496;
pktTotal        = 8 + pktLen;          % 1504 bits per TS packet
numPktsPerFrame = wavegen.MinNumPackets(1);
dataSize        = numPktsPerFrame * pktTotal;

%% =========================================================
%  PHY FRAME GEOMETRY (fixed by MODCOD)
%% =========================================================
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

% Reference pilots (scrambled, unmodulated)
plScrambIntSeq = satcom.internal.dvbs.plScramblingIntegerSequence(0);
cMap           = [1; 1j; -1; -1j];
cSeq           = cMap(plScrambIntSeq + 1);
[~, pilotInd_raw] = satcom.internal.dvbs.pilotBlock(numPilotBlks);
pilotInd   = pilotInd_raw + slotLen;
refPilots  = (1+1j)/sqrt(2) .* cSeq(pilotInd_raw);

filterDelay = wavegen.FilterSpanInSymbols / 2;
searchRange = 5;

%% =========================================================
%  MAIN Es/No SWEEP
%% =========================================================
for i = 1:numEsNo

    bitsErr       = 0;
    pktsErr       = 0;
    pktsRec       = 0;
    nmseAcc       = 0;    % accumulate NMSE across frames
    nmseCount     = 0;
    numFramesLost = 0;
    totalFrames   = 0;
    batchSize     = 20;
    batchStart    = 1;

    % Pre-generate bit pool for all maxFrames at this Es/No point
    numPkts   = numPktsPerFrame * maxFrames;
    pn = comm.PNSequence('Polynomial','x9+x5+1', ...
        'InitialConditions', [zeros(1,8) 1], ...
        'VariableSizeOutput', true, ...
        'MaximumOutputSize',  [pktLen * numPkts, 1]);
    reset(pn);
    rawpkts   = pn(pktLen * numPkts);
    txRawPkts = reshape(rawpkts, pktLen, numPkts);
    txPkts    = [repmat(syncBits, 1, numPkts); txRawPkts];
    inputBits = txPkts(:);

    % snr_samp: per-sample SNR accounting for oversampling (SPS=2)
    % This is used to set noise power relative to received signal power.
    snr_samp = esno(i) - 10*log10(wavegen.SamplesPerSymbol);

    % ── Adaptive batch loop ────────────────────────────────────────────────
    while bitsErr < minBitErrors && totalFrames < maxFrames

        framesThisBatch = min(batchSize, maxFrames - totalFrames);
        pktsThisBatch   = numPktsPerFrame * framesThisBatch;
        pktOffset       = (batchStart - 1) * numPktsPerFrame;

        batchBits = inputBits((pktOffset*pktTotal + 1) : ...
                              (pktOffset + pktsThisBatch)*pktTotal);
        batchData = {batchBits};

        % ── Transmitter ───────────────────────────────────────────────────
        release(wavegen);
        wavegen.DFL = getdfl(wavegen.MODCOD, wavegen.FECFrame);
        txWaveform  = [wavegen(batchData); flushFilter(wavegen)];

        % ── Rain fading channel ───────────────────────────────────────────
        % Apply the pre-computed channel coefficient h_k to all samples.
        % A new random phase is drawn per batch to model slow fading
        % varying across frames (one realisation per batch ≈ coherence time).
        phi_batch  = 2*pi*rand(1);
        h_k_batch  = (1/Ap^0.5) * exp(-1j*phi_batch) * sqrt(b_gain) * sqrt(b_max);
        attWaveform = h_k_batch * txWaveform;

        % ── AWGN — noise power referenced to POST-FADE signal ─────────────
        % rxPower is the actual signal power at the receiver input.
        % This ensures Es/N0 is correct at the decoder, regardless of
        % how much h_k has attenuated the signal.
        rxPower  = mean(abs(attWaveform).^2);
        noiseVar = rxPower / (10^(snr_samp/10));
        noise    = sqrt(noiseVar/2) * ...
                   (randn(size(attWaveform)) + 1j*randn(size(attWaveform)));
        rxWaveform = attWaveform + noise;

        % ── Matched filter (RRC receive filter) ───────────────────────────
        rxFilt = comm.RaisedCosineReceiveFilter( ...
            'RolloffFactor',         wavegen.RolloffFactor, ...
            'FilterSpanInSymbols',   wavegen.FilterSpanInSymbols, ...
            'InputSamplesPerSymbol', wavegen.SamplesPerSymbol, ...
            'DecimationFactor',      wavegen.SamplesPerSymbol);
        rxSym_raw = rxFilt(rxWaveform);
        release(rxFilt);

        % Trim filter group-delay transient (filterDelay symbols each side)
        rxSymbols = rxSym_raw(filterDelay + 1 : end - filterDelay);
        % DO NOT normalise here — LS estimation needs the raw channel-scaled
        % symbols so that H_ls captures the true h_k magnitude.

        % ── Per-frame processing ──────────────────────────────────────────
        for stIdx = 1:framesThisBatch

            if bitsErr >= minBitErrors || totalFrames >= maxFrames
                break;
            end

            nominalStart = (stIdx-1)*plFrameSize + 1;
            nominalEnd   = nominalStart + plFrameSize - 1;
            if nominalEnd + searchRange > length(rxSymbols)
                break;
            end

            % ── Coarse frame sync (pilot correlation on raw symbols) ──────
            % Performed BEFORE equalization — correlation metric is
            % phase-invariant: abs(mean(y_p .* p*)) does not require
            % knowing the channel phase.
            bestOffset = 0;
            bestMetric = -Inf;
            for offset = -searchRange:searchRange
                cStart = nominalStart + offset;
                if cStart < 1 || (cStart + plFrameSize - 1) > length(rxSymbols)
                    continue;
                end
                cand   = rxSymbols(cStart : cStart + plFrameSize - 1);
                metric = abs(mean(cand(pilotInd) .* conj(refPilots)));
                if metric > bestMetric
                    bestMetric = metric;
                    bestOffset = offset;
                end
            end

            frameStart = nominalStart + bestOffset;
            rxFrame    = rxSymbols(frameStart : frameStart + plFrameSize - 1);
            totalFrames = totalFrames + 1;

            % ── LS Channel Estimation (Eq. 12–16) ─────────────────────────
            % Estimate h_k from all pilot positions in this frame:
            %   ĥ_LS = (1/N_p) Σ y_p[n] / p[n]
            % Average within each pilot block first, then across blocks
            % to reduce noise.
            rxPilSym   = rxFrame(pilotInd);              % [pilotSymLen×numBlks, 1]
            LS_raw     = rxPilSym ./ refPilots;           % element-wise division
            LS_mat     = reshape(LS_raw, pilotSymLen, numPilotBlks);
            LS_perBlk  = mean(LS_mat, 1);                 % [1 × numPilotBlks]
            H_ls       = mean(LS_perBlk);                 % scalar estimate

            % NMSE = |h_k_batch - H_ls|² / |h_k_batch|²
            nmseFrame = abs(h_k_batch - H_ls)^2 / abs(h_k_batch)^2;
            nmseAcc   = nmseAcc + nmseFrame;
            nmseCount = nmseCount + 1;

            % ── Equalization ──────────────────────────────────────────────
            rxFrameEq = rxFrame ./ H_ls;

            % ── PL Header Recovery ────────────────────────────────────────
            phyParams = dvbsPLHeaderRecover(rxFrameEq(1:slotLen), ...
                                            Mode="DVB-S2/S2X regular");
            M = phyParams.ModulationOrder;
            if M == 0; R = []; else; R = eval(phyParams.LDPCCodeIdentifier); end

            if M ~= modOrder || R ~= codeRate || ...
                    phyParams.FECFrameLength ~= cwLen || ~phyParams.HasPilots
                numFramesLost = numFramesLost + 1;
                continue;
            end

            % ── Noise Variance Estimate (pilot-aided, post-equalization) ──
            nVar = DVBS2NoiseVarEstimate(rxFrameEq, pilotInd, refPilots, normFlag);

            % ── FEC Decoding ──────────────────────────────────────────────
            [decBitsTemp, isFrameLost, pktCRC] = dvbs2BitRecover( ...
                rxFrameEq, nVar, normFlag);
            decBits = decBitsTemp{:};

            if ~isFrameLost && length(decBits) ~= dataSize
                isFrameLost = true;
            end
            numFramesLost = numFramesLost + isFrameLost;

            if ~isFrameLost
                % PER
                pktsErr = pktsErr + numel(pktCRC{:}) - sum(pktCRC{:});
                pktsRec = pktsRec + numel(pktCRC{:});

                % BER
                globalPktOffset = pktOffset + (stIdx-1)*numPktsPerFrame;
                bitInd = globalPktOffset*pktTotal + (1:dataSize);
                if bitInd(end) <= length(inputBits)
                    bitsErr = bitsErr + sum(inputBits(bitInd) ~= decBits);
                end
            end
        end % per-frame loop

        batchStart = batchStart + framesThisBatch;
        if batchStart * numPktsPerFrame > numPkts
            break;
        end

    end % adaptive batch loop

    % ── Aggregate results ─────────────────────────────────────────────────
    if pktsRec > 0
        berResults(i)  = bitsErr / (pktsRec * pktTotal);
        perResults(i)  = pktsErr / pktsRec;
    else
        berResults(i)  = 1.0;
        perResults(i)  = 1.0;
    end

    if nmseCount > 0
        nmseResults(i) = nmseAcc / nmseCount;
    end

    confFlag = '';
    if bitsErr < minBitErrors
        confFlag = '  [<200 errors — low confidence]';
    end

    fprintf(['Es/No = %5.2f dB | ChanPwr = %.2f dB | ' ...
             'Frames %d/%d lost | Errors %d | ' ...
             'BER = %.3e | PER = %.3e | NMSE = %.3e%s\n'], ...
        esno(i), chanPow_dB, numFramesLost, totalFrames, bitsErr, ...
        berResults(i), perResults(i), nmseResults(i), confFlag);

    clearvars -except esno numEsNo berResults perResults nmseResults i ...
        wavegen modOrder codeRate cwLen dataLen slotLen pilotLen plFrameSize ...
        numPilotBlks pilotInd pilotInd_raw refPilots cSeq ...
        numPktsPerFrame dataSize pktTotal pktLen syncBits normFlag ...
        filterDelay searchRange minBitErrors maxFrames ...
        chanParams sysP pRain Ap_dB Ap b_max b_gain chanPow_dB;
end

%% =========================================================
%  SAVE
%% =========================================================
save('dvbs2_ls_rainfade_results.mat', ...
    'esno', 'berResults', 'perResults', 'nmseResults', ...
    'Ap_dB', 'chanPow_dB');
fprintf('\nSaved → dvbs2_ls_rainfade_results.mat\n');

%% =========================================================
%  PLOT
%% =========================================================
validBER  = berResults  > 0 & ~isnan(berResults);
validPER  = perResults  > 0 & ~isnan(perResults);
validNMSE = nmseResults > 0 & ~isnan(nmseResults);

figure('Name','DVB-S2 LS Rain Fading — BER & PER');
semilogy(esno(validBER), berResults(validBER),  'b-o', ...
    'LineWidth',1.5,'MarkerFaceColor','b'); hold on;
semilogy(esno(validPER), perResults(validPER),  'r-s', ...
    'LineWidth',1.5,'MarkerFaceColor','r');
grid on;
legend('BER (LS est.)','PER (LS est.)','Location','southwest');
xlabel('E_s/N_0 at receiver input (dB)');
ylabel('Error rate');
title(sprintf(['DVB-S2 BER & PER — LS Rain Fading Channel\n' ...
    'MODCOD %d | Tropical (Penang) | A_p = %.2f dB'], ...
    wavegen.MODCOD, Ap_dB));
ylim([1e-6 1]);

figure('Name','DVB-S2 LS Rain Fading — NMSE');
semilogy(esno(validNMSE), nmseResults(validNMSE), 'g-^', ...
    'LineWidth',1.5,'MarkerFaceColor','g');
grid on;
xlabel('E_s/N_0 at receiver input (dB)');
ylabel('NMSE');
title(sprintf(['LS Channel Estimation NMSE\n' ...
    'MODCOD %d | Tropical (Penang) | A_p = %.2f dB'], ...
    wavegen.MODCOD, Ap_dB));

%% =========================================================
%  HELPER FUNCTIONS
%% =========================================================

function dfl = getdfl(modCod, fecFrame)
%GETDFL  Data field length lookup table.
    if strcmp(fecFrame,'normal')
        lut = [16008 21408 25728 32208 38688 43040 48408 51648 53840 57472 ...
               58192 38688 43040 48408 53840 57472 58192 43040 48408 51648 ...
               53840 57472 58192 48408 51648 53840 57472 58192] - 80;
    else
        lut = [3072 5232 6312 7032 9552 10632 11712 12432 13152 14232 0 ...
               9552 10632 11712 13152 14232 0 10632 11712 12432 13152 14232 ...
               0 11712 12432 13152 14232 0] - 80;
    end
    dfl = lut(modCod);
end

% -----------------------------------------------------------------

function nVarEst = DVBS2NoiseVarEstimate(rxData, pilotInd, refPilots, normFlag)
%DVBS2NOISEVARESTIMATE  Pilot-aided noise variance estimator.
%
%  Subtracts the coherent signal power from total pilot power to
%  isolate the noise.  normFlag passed from ber_analysis convention
%  (applied to rxData before estimation when true).
    if normFlag
        rxData = rxData / sqrt(mean(abs(rxData).^2));
    end
    rxPilots = rxData(pilotInd, :);
    pSigN    = mean(abs(rxPilots).^2);                    % E[|y|²]
    pSig     = abs(mean(rxPilots .* conj(refPilots))).^2; % |E[y·p*]|²
    nVarEst  = max(pSigN - pSig, 0);                      % clamp to ≥0
end

% -----------------------------------------------------------------

function Ap_dB = calcRainAttenuationP618(chanParams, p)
%CALCRAINATTENUATIONP618  ITU-R P.618-13 rain attenuation.
    f_GHz = chanParams.Frequency_GHz;
    theta  = chanParams.ElevationAngle;
    lat    = chanParams.Latitude;
    hs_km  = chanParams.Altitude_m / 1000;
    R001   = chanParams.R001_mmph;
    h0_km  = chanParams.IsothermHeight_km;
    hR_km  = h0_km + 0.36;
    Re_km  = 8500;

    if R001 <= 0 || hR_km <= hs_km
        Ap_dB = 0; return;
    end

    if theta >= 5
        Ls = (hR_km - hs_km) / sind(theta);
    else
        Ls = 2*(hR_km-hs_km) / ...
             (sqrt(sind(theta)^2 + 2*(hR_km-hs_km)/Re_km) + sind(theta));
    end
    LG = Ls * cosd(theta);

    [k_c, alpha_c] = ituP838_k_alpha(f_GHz, theta, chanParams.Polarization);
    gamma_R = k_c * (R001^alpha_c);

    r001  = 1 / (1 + 0.78*sqrt(LG*gamma_R/f_GHz) - 0.38*(1-exp(-2*LG)));
    zeta  = atand((hR_km-hs_km) / (LG*r001));
    if zeta > theta
        LR = (LG*r001) / cosd(theta);
    else
        LR = (hR_km-hs_km) / sind(theta);
    end

    chi  = max(36 - abs(lat), 0);
    v001 = 1 / (1 + sqrt(sind(theta)) * ...
           (31*(1-exp(-theta/(1+chi)))*sqrt(LR*gamma_R)/f_GHz^2 - 0.45));

    LE   = LR * v001;
    A001 = gamma_R * LE;

    if abs(p - 0.01) < 1e-12
        Ap_dB = A001; return;
    end

    if p >= 1 || abs(lat) >= 36
        beta = 0;
    elseif theta >= 25
        beta = -0.005*(abs(lat)-36);
    else
        beta = -0.005*(abs(lat)-36) + 1.8 - 4.25*sind(theta);
    end

    if A001 <= 0
        Ap_dB = 0; return;
    end

    exp_val = -(0.655 + 0.033*log(p) - 0.045*log(A001) ...
               - beta*(1-p)*sind(theta));
    Ap_dB = A001 * (p/0.01)^exp_val;
end

% -----------------------------------------------------------------

function [k, alpha] = ituP838_k_alpha(f_GHz, theta_deg, polarization)
%ITUP838_K_ALPHA  Specific attenuation coefficients (ITU-R P.838-3).
    x = log10(f_GHz);

    akH=[-5.33980,-0.35351,-0.23789,-0.94158]; bkH=[-0.10008,1.26970,0.86036,0.64552];
    ckH=[1.13098,0.45400,0.15354,0.16817]; mkH=-0.18961; ckH0=0.71147;
    logkH = sum(akH.*exp(-((x-bkH)./ckH).^2)) + mkH*x + ckH0;
    kH    = 10^logkH;

    akV=[-3.80595,-3.44965,-0.39902,0.50167]; bkV=[0.56934,-0.22911,0.73042,1.07319];
    ckV=[0.81061,0.51059,0.11899,0.27195]; mkV=-0.16398; ckV0=0.63297;
    logkV = sum(akV.*exp(-((x-bkV)./ckV).^2)) + mkV*x + ckV0;
    kV    = 10^logkV;

    aaH=[-0.14318,0.29591,0.32177,-5.37610,16.1721];
    baH=[1.82442,0.77564,0.63773,-0.96230,-3.29980];
    caH=[-0.55187,0.19822,0.13164,1.47828,3.43990];
    maH=0.67849; caH0=-1.95537;
    alphaH = sum(aaH.*exp(-((x-baH)./caH).^2)) + maH*x + caH0;

    aaV=[-0.07771,0.56727,-0.20238,-48.2991,48.5833];
    baV=[2.33840,0.95545,1.14520,0.791669,0.791459];
    caV=[-0.76284,0.54039,0.26809,0.116226,0.116479];
    maV=-0.053739; caV0=0.83433;
    alphaV = sum(aaV.*exp(-((x-baV)./caV).^2)) + maV*x + caV0;

    pol = upper(string(polarization));
    if pol == "H"; tau = 0; elseif pol == "V"; tau = 90; else; tau = 45; end

    k     = (kH + kV + (kH-kV)*cosd(theta_deg)^2*cosd(2*tau)) / 2;
    alpha = (kH*alphaH + kV*alphaV + ...
             (kH*alphaH-kV*alphaV)*cosd(theta_deg)^2*cosd(2*tau)) / (2*k);
end
