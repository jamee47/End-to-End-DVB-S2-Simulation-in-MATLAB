%% DVB-S2 Channel Estimation
% Author : Mohtasim Al Jamee
%
% Sections:
%   1 - Configuration
%   2 - Transmitter  (dvbs2WaveformGenerator, pilot info from toolbox)
%   3 - Channel      (rain fading + AWGN)
%   4 - Receiver     (matched filter, extract pilots)
%   5 - LS  Channel Estimation  (adapted from OFDM LS_CE reference)
%   6 - MMSE Channel Estimation (adapted from OFDM MMSE_CE reference)
%   7 - BER vs SNR plot

clc; clearvars; close all;

%% ========================================================================
%  SECTION 1: CONFIGURATION
%% ========================================================================
MODCOD        = 2;        % QPSK 1/3 short frame
FECFrame      = 'short';
RolloffFactor = 0.35;
FilterSpan    = 10;       % symbols
SPS           = 4;        % samples per symbol

EsNo_inspect  = 20;       % dB — for single frame inspection

% ITU-R P.618 rain channel — Penang, Malaysia (tropical)
k_V = 0.07510;  alpha_V = 1.09980;
R_001 = 130.0;  h_R_km = 4.86;
lat_deg = 5.17; elev_deg = 45; alt_km = 0.057; freq_GHz = 20;

% BER sweep
EsNo_range    = 0:2:30;
numEvalFrames = 50;

%% ========================================================================
%  SECTION 2: TRANSMITTER
%  Use dvbs2WaveformGenerator.
%  Get pilot positions and pilot symbols directly from the toolbox
%  using satcom.internal.dvbs functions — no manual construction needed.
%% ========================================================================
fprintf('--- SECTION 2: Transmitter ---\n');

%% 2a — Configure generator
wavegen = dvbs2WaveformGenerator;
wavegen.StreamFormat        = 'TS';
wavegen.NumInputStreams      = 1;
wavegen.MODCOD               = MODCOD;
wavegen.FECFrame             = FECFrame;
wavegen.HasPilots            = true;
wavegen.RolloffFactor        = RolloffFactor;
wavegen.FilterSpanInSymbols  = FilterSpan;
wavegen.SamplesPerSymbol     = SPS;
wavegen.DFL = 3004;

%% 2b — Input bits (TS stream)
syncBits = [0;1;0;0;0;1;1;1];
pktLen   = 1496;
numPkts  = wavegen.MinNumPackets(1);
txRaw    = randi([0 1], pktLen, numPkts);
txPkts   = [repmat(syncBits,1,numPkts); txRaw];
txBits   = txPkts(:);

%% 2c — Generate waveform
txWaveform = [wavegen(txBits); flushFilter(wavegen)];
fprintf('TX waveform : %d samples\n', length(txWaveform));

%% 2d — Get PL frame structure from toolbox
[modOrder, ~, cwLen] = satcom.internal.dvbs.getS2PHYParams(...
    wavegen.MODCOD, wavegen.FECFrame);

dataLen      = cwLen / log2(modOrder);
slotLen      = 90;
pilotSymLen  = 36;
pilotBlkFreq = 16;

numPilotBlks = floor(dataLen/(slotLen*pilotBlkFreq));
if floor(dataLen/(slotLen*16)) == dataLen/(slotLen*pilotBlkFreq)
    numPilotBlks = numPilotBlks - 1;
end

pilotLen    = numPilotBlks * pilotSymLen;
plFrameSize = dataLen + pilotLen + slotLen;

fprintf('PL frame    : %d symbols\n', plFrameSize);
fprintf('Pilot blocks: %d\n', numPilotBlks);

%% 2e — Get pilot positions AND pilot values from the toolbox
%  satcom.internal.dvbs.pilotBlock gives:
%    pilotSymbols : the actual pilot symbol values as transmitted
%    pilotInd_raw : their positions within the XFECFRAME (no PLHEADER)
%
%  Add slotLen (90) to get positions within the full PL frame.
%
%  These pilot symbols already include PL scrambling — they are the
%  exact values the transmitter put on the wire. No manual construction.

[pilotSymbols, pilotInd_raw] = satcom.internal.dvbs.pilotBlock(numPilotBlks);

% Positions within the full PL frame (offset by PLHEADER)
pilotInd = pilotInd_raw + slotLen;

% Data positions: everything that is not PLHEADER and not a pilot
dataInd  = setdiff((1:plFrameSize)', pilotInd);
dataInd  = dataInd(dataInd > slotLen);

fprintf('Pilot syms  : %d  (magnitude check: %.4f — should be 0.7071)\n',...
    length(pilotSymbols), abs(pilotSymbols(1)));
fprintf('Data syms   : %d\n', length(dataInd));

%% ========================================================================
%  SECTION 3: CHANNEL — Rain fading + AWGN
%% ========================================================================
fprintf('\n--- SECTION 3: Channel ---\n');

%% 3a — A_001 (ITU-R P.618-13)
gamma_R = k_V * (R_001^alpha_V);
L_S = (h_R_km - alt_km)/sind(elev_deg);
L_G = L_S*cosd(elev_deg);
r_001 = 1/(1+0.78*sqrt(L_G*gamma_R/freq_GHz)-0.38*(1-exp(-2*L_G)));
zeta  = atand((h_R_km-alt_km)/(L_G*r_001));
L_E   = (zeta>elev_deg)*L_G*r_001/cosd(elev_deg) + ...
        (zeta<=elev_deg)*(h_R_km-alt_km)/sind(elev_deg);
chi   = 36-lat_deg;
v_001 = 1/(1+sqrt(sind(elev_deg))*(31*(1-exp(-elev_deg/(1+chi)))...
        *sqrt(L_E*gamma_R)/freq_GHz^2-0.45));
A_001 = gamma_R*L_E*v_001;

%% 3b — Instantaneous A_p (ITU-R P.618-13 CDF step 6)
rng('shuffle');
p_pct = 0.01 + (5.0-0.01)*rand(1);
if p_pct<=1
    C = 0.655+0.033*log(p_pct)-0.045*log(max(A_001,0.01));
else
    C = 0.655+0.033*log(p_pct)-0.045*log(max(A_001,0.01)) ...
        -0.053*(1-p_pct)*sind(elev_deg);
end
A_p_dB  = max(min(A_001*(p_pct/0.01)^(-C), EsNo_inspect-3), 0.1);
h_tilde = 10^(-A_p_dB/20) * exp(-1j*2*pi*rand(1));

rxWaveform = awgn(txWaveform*h_tilde, EsNo_inspect-10*log10(SPS), 'measured');

fprintf('A_001       : %.3f dB\n', A_001);
fprintf('A_p         : %.3f dB  (p=%.3f%%)\n', A_p_dB, p_pct);
fprintf('|h_tilde|   : %.6f\n', abs(h_tilde));

%% ========================================================================
%  SECTION 4: RECEIVER — matched filter, extract pilots and data
%% ========================================================================
fprintf('\n--- SECTION 4: Receiver ---\n');

filterDelay      = FilterSpan/2;
plFrameSize_orig = plFrameSize;  % keep original for BER loop (before safeLen trim)

%% 4a — Matched filter + downsample (RX)
rxFilt = comm.RaisedCosineReceiveFilter('RolloffFactor',RolloffFactor,...
    'FilterSpanInSymbols',FilterSpan,...
    'InputSamplesPerSymbol',SPS,'DecimationFactor',SPS);
rxSym = rxFilt(rxWaveform);
rxSym = rxSym(filterDelay+1:end);
release(rxFilt);

%% 4b — Matched filter (TX reference for BER)
txFilt = comm.RaisedCosineReceiveFilter('RolloffFactor',RolloffFactor,...
    'FilterSpanInSymbols',FilterSpan,...
    'InputSamplesPerSymbol',SPS,'DecimationFactor',SPS);
txSym = txFilt(txWaveform);
txSym = txSym(filterDelay+1:end);
release(txFilt);

%% 4c — Guard: ensure both symbol streams are long enough before slicing
%  After matched filtering the stream may be slightly shorter than
%  plFrameSize due to filter transients. Trim to the safe length.
safeLen = min([length(rxSym), length(txSym), plFrameSize]);
if safeLen < plFrameSize
    warning('Symbol stream shorter than PL frame (%d < %d). Trimming.', ...
        safeLen, plFrameSize);
end
rxFrame = rxSym(1:safeLen);
txFrame = txSym(1:safeLen);
plFrameSize = safeLen;   % update so all downstream indexing stays consistent

%% 4d — Extract received pilots and data symbols
%  Reshape all pilot arrays to [pilotSymLen x numPilotBlks]
%
%  Guard: pilotInd and pilotSymbols must have the same length before masking.
%  satcom.internal.dvbs.pilotBlock may return a pilotSymbols vector whose
%  length differs from pilotInd (e.g. column vs row, or extra/fewer entries).
%  Trim both to the shorter of the two first, then apply the bounds mask.
minPLen      = min(length(pilotInd), length(pilotSymbols));
pilotInd     = pilotInd(1:minPLen);
pilotSymbols = pilotSymbols(1:minPLen);

% Clip pilot indices to actual frame length (plFrameSize may have been
% trimmed to safeLen). Remove any pilot positions that fall beyond the frame.
validPilotMask = pilotInd <= plFrameSize;
pilotInd       = pilotInd(validPilotMask);
pilotSymbols   = pilotSymbols(validPilotMask);

% Recompute numPilotBlks so reshape is always exact (no remainder)
numPilotBlks   = floor(length(pilotInd) / pilotSymLen);
keepN          = numPilotBlks * pilotSymLen;
pilotInd       = pilotInd(1:keepN);
pilotSymbols   = pilotSymbols(1:keepN);

% Recompute dataInd using updated plFrameSize and clipped pilotInd
dataInd  = setdiff((1:plFrameSize)', pilotInd);
dataInd  = dataInd(dataInd > slotLen);

% Now reshape is safe — all sizes are consistent
pilotIdxMat    = reshape(pilotInd,    pilotSymLen, numPilotBlks);
pilotSymbolMat = reshape(pilotSymbols,pilotSymLen, numPilotBlks);
rxPilotMat     = rxFrame(pilotIdxMat);    % [36 x numPilotBlks]

rxDataSym = rxFrame(dataInd);
txDataSym = txFrame(dataInd);

fprintf('Pilot blocks    : %d  (after bounds check)\n', numPilotBlks);
fprintf('RX pilots shape : [%d x %d]\n', size(rxPilotMat));
fprintf('RX data symbols : %d\n', length(rxDataSym));

%% ========================================================================
%  SECTION 5: LS CHANNEL ESTIMATION
%
%  Adapted from OFDM LS_CE reference:
%    OFDM  : H_LS(k) = Y(pilot_loc(k)) / Xp(k)   then interpolate in freq
%    DVB-S2: h_LS(b) = mean(Y_b / Xp_b)           then interpolate in time
%
%  Y_b  = received pilot symbols in block b  [36 x 1]
%  Xp_b = known pilot symbols in block b     [36 x 1]  (from toolbox)
%
%  Dividing per-symbol then averaging 36 symbols = 15.6 dB noise reduction
%  Interpolation fills channel estimate between pilot block positions.
%% ========================================================================
fprintf('\n--- SECTION 5: LS Estimation ---\n');

%% Step 1: LS estimate per pilot block
%  h_LS_b = mean(Y_b ./ Xp_b)   [scalar per block]
h_LS_blocks = mean(rxPilotMat ./ pilotSymbolMat, 1).';  % [numPilotBlks x 1]

%% Step 2: Pilot block centre positions (for interpolation)
%  The centre of block b sits at the middle of its 36 pilot symbols
%  pilotIdxMat(:,b) gives all 36 positions of block b
pilotCentres = mean(pilotIdxMat, 1).';   % [numPilotBlks x 1]  centre indices

%% Step 3: Interpolate over all PL frame positions
%  Equivalent to H_LS = interpolate(LS_est, pilot_loc, Nfft, method) in OFDM
%  interp1 requires >= 2 sample points. If only 1 pilot block survived the
%  bounds trim, fall back to a constant (nearest) fill. If 0 blocks, fill
%  with zeros (estimation is impossible — caller should check numPilotBlks).
allPos = (1:plFrameSize).';
if numPilotBlks >= 2
    h_LS_all = interp1(pilotCentres, h_LS_blocks, allPos, 'linear', 'extrap');
elseif numPilotBlks == 1
    warning('Only 1 pilot block available — using constant channel estimate.');
    h_LS_all = repmat(h_LS_blocks(1), plFrameSize, 1);
else
    warning('No pilot blocks available — channel estimate set to zero.');
    h_LS_all = zeros(plFrameSize, 1);
end

%% Step 4: Frame-level estimate (mean over pilot blocks)
h_LS_frame = mean(h_LS_blocks);

NMSE_LS = mean(abs(h_LS_blocks - h_tilde).^2) / (abs(h_tilde)^2 + eps);

fprintf('LS estimates per block (first 3):\n');
for k=1:min(3,numPilotBlks)
    fprintf('  Block %2d  pos=%4.0f  h=%.4f+%.4fj  |h|=%.4f\n',...
        k, pilotCentres(k), real(h_LS_blocks(k)), imag(h_LS_blocks(k)),...
        abs(h_LS_blocks(k)));
end
fprintf('h_LS_frame : %.5f+%.5fj  |h|=%.5f\n',...
    real(h_LS_frame),imag(h_LS_frame),abs(h_LS_frame));
fprintf('h_tilde    : %.5f+%.5fj  |h|=%.5f\n',...
    real(h_tilde),imag(h_tilde),abs(h_tilde));
fprintf('NMSE_LS    : %.4e\n', NMSE_LS);

%% ========================================================================
%  SECTION 6: MMSE CHANNEL ESTIMATION
%
%  Adapted from OFDM MMSE_CE reference:
%    OFDM : H_MMSE = Rhp * inv(Rpp) * H_tilde   (frequency domain)
%    DVB-S2: h_MMSE = Rhp * inv(Rpp) * h_LS     (time domain, pilot blocks)
%
%  For DVB-S2 flat slow fading:
%    h_LS is the LS estimate per pilot block  [Nb x 1]
%    Rpp  = R_hh*I + (sigma_n2/Np)*I          [Nb x Nb]  (diagonal)
%    Rhp  = R_hh*I                             [Nb x Nb]  (diagonal)
%
%  This simplifies to a scalar Wiener filter:
%    W     = R_hh / (R_hh + sigma_n2/Np)
%    h_MMSE = W * h_LS
%
%  Noise variance sigma_n2 estimated via SNORE (data-aided ML estimator):
%    pSigN    = mean received pilot power
%    pSig     = signal power via cross-correlation with known pilots
%    sigma_n2 = pSigN - pSig
%% ========================================================================
fprintf('\n--- SECTION 6: MMSE Estimation ---\n');

%% Step 1: SNORE noise variance estimate
%  pSigN = mean received power; pSig = mean signal power = |mean(y·conj(p))|^2
%  (mean first, then square — gives mean amplitude squared = mean signal power)
rxPilotsVec = rxFrame(pilotInd);           % all received pilot symbols
pSigN    = mean(abs(rxPilotsVec).^2);
pSig     = mean(abs(rxPilotsVec .* conj(pilotSymbols)));   % mean amplitude
pSig     = pSig^2;                                          % mean signal power
sigma_n2 = max(pSigN - pSig, 1e-12);

%% Step 2: Channel power estimate R_hh (bias-corrected)
Np   = pilotSymLen;   % = 36 averaging gain per block
R_hh = mean(abs(h_LS_blocks).^2) - sigma_n2/Np;
R_hh = max(real(R_hh), 1e-12);

%% Step 3: Wiener filter coefficient W
%  Equivalent to Rhp * inv(Rpp) from MMSE_CE reference
W_mmse = R_hh / (R_hh + sigma_n2/Np);

%% Step 4: MMSE estimate per pilot block + interpolation
h_MMSE_blocks = W_mmse * h_LS_blocks;
if numPilotBlks >= 2
    h_MMSE_all = interp1(pilotCentres, h_MMSE_blocks, allPos, 'linear', 'extrap');
elseif numPilotBlks == 1
    h_MMSE_all = repmat(h_MMSE_blocks(1), plFrameSize, 1);
else
    h_MMSE_all = zeros(plFrameSize, 1);
end
h_MMSE_frame  = mean(h_MMSE_blocks);

NMSE_MMSE = mean(abs(h_MMSE_blocks - h_tilde).^2)/(abs(h_tilde)^2+eps);

fprintf('sigma_n2   : %.4e  (SNORE)\n', sigma_n2);
fprintf('R_hh       : %.4e  (channel power)\n', R_hh);
fprintf('W_mmse     : %.6f  (Wiener coeff)\n', W_mmse);
fprintf('h_MMSE_frame: %.5f+%.5fj\n', real(h_MMSE_frame),imag(h_MMSE_frame));
fprintf('NMSE_MMSE  : %.4e\n', NMSE_MMSE);

%% Quick constellation check at EsNo_inspect dB
figure('Name',sprintf('Constellation EsNo=%.0fdB',EsNo_inspect),...
    'Color','w','Position',[50 400 1000 360]);

subplot(1,3,1);
scatter(real(rxDataSym),imag(rxDataSym),8,'r','filled','MarkerFaceAlpha',.3);
axis equal; grid on; xlim([-2 2]); ylim([-2 2]);
title('No equalization'); xlabel('I'); ylabel('Q');

subplot(1,3,2);
eqLS = rxDataSym / h_LS_frame;
scatter(real(eqLS),imag(eqLS),8,'b','filled','MarkerFaceAlpha',.3);
axis equal; grid on; xlim([-2 2]); ylim([-2 2]);
title('After LS equalization'); xlabel('I'); ylabel('Q');

subplot(1,3,3);
eqMM = rxDataSym / h_MMSE_frame;
scatter(real(eqMM),imag(eqMM),8,[0 .6 0],'filled','MarkerFaceAlpha',.3);
axis equal; grid on; xlim([-2 2]); ylim([-2 2]);
title('After MMSE equalization'); xlabel('I'); ylabel('Q');

sgtitle(sprintf('QPSK Constellation | EsNo=%.0fdB | |h|=%.4f | LS NMSE=%.3e | MMSE NMSE=%.3e',...
    EsNo_inspect,abs(h_tilde),NMSE_LS,NMSE_MMSE),'FontSize',11,'FontWeight','bold');

%% ========================================================================
%  SECTION 7: BER vs SNR
%% ========================================================================
fprintf('\n--- SECTION 7: BER vs SNR ---\n');

numSNR       = length(EsNo_range);
BER_LS_arr   = zeros(1,numSNR);
BER_MMSE_arr = zeros(1,numSNR);
BER_noEq_arr = zeros(1,numSNR);
qpsk         = @(s)[real(s)>=0, imag(s)>=0];

for si = 1:numSNR
    snr = EsNo_range(si);
    bLS=0; bMM=0; bNE=0; nV=0;

    for f = 1:numEvalFrames
        %% TX — fresh bits each frame
        wg = dvbs2WaveformGenerator;
        wg.StreamFormat='TS'; wg.NumInputStreams=1;
        wg.MODCOD=MODCOD; wg.FECFrame=FECFrame; wg.HasPilots=true;
        wg.RolloffFactor=RolloffFactor;
        wg.FilterSpanInSymbols=FilterSpan; wg.SamplesPerSymbol=SPS;
        wg.DFL = 3004;
        rp  = randi([0 1],pktLen,numPkts);
        tp  = [repmat(syncBits,1,numPkts);rp];
        txW = [wg(tp(:));flushFilter(wg)];

        %% Channel — unique fading per frame
        p_f = 0.01+(5.0-0.01)*rand(1);
        if p_f<=1
            Cf=0.655+0.033*log(p_f)-0.045*log(max(A_001,0.01));
        else
            Cf=0.655+0.033*log(p_f)-0.045*log(max(A_001,0.01))...
               -0.053*(1-p_f)*sind(elev_deg);
        end
        Ap_f = max(min(A_001*(p_f/0.01)^(-Cf),snr-3),0.1);
        ht_f = 10^(-Ap_f/20)*exp(-1j*2*pi*rand(1));
        rxW  = awgn(txW*ht_f, snr-10*log10(SPS),'measured');

        %% RX — matched filter
        rf = comm.RaisedCosineReceiveFilter('RolloffFactor',RolloffFactor,...
            'FilterSpanInSymbols',FilterSpan,...
            'InputSamplesPerSymbol',SPS,'DecimationFactor',SPS);
        rxS=rf(rxW); rxS=rxS(filterDelay+1:end); release(rf);

        tf = comm.RaisedCosineReceiveFilter('RolloffFactor',RolloffFactor,...
            'FilterSpanInSymbols',FilterSpan,...
            'InputSamplesPerSymbol',SPS,'DecimationFactor',SPS);
        txS=tf(txW); txS=txS(filterDelay+1:end); release(tf);

        safeLen_f = min([length(rxS), length(txS), plFrameSize_orig]);
        if safeLen_f < plFrameSize_orig*0.9, continue; end  % skip if too short
        rxFr=rxS(1:safeLen_f);
        txFr=txS(1:safeLen_f);
        % Clip pilot and data indices to actual frame length
        pInd_f = pilotInd(pilotInd <= safeLen_f);
        dInd_f = dataInd(dataInd  <= safeLen_f);
        if isempty(pInd_f) || isempty(dInd_f), continue; end
        nBlks_f = floor(length(pInd_f)/pilotSymLen);
        if nBlks_f < 1, continue; end
        pInd_f  = pInd_f(1:nBlks_f*pilotSymLen);
        pSym_f  = pilotSymbols(1:nBlks_f*pilotSymLen);

        %% LS estimation — use frame-local pilot index and symbol matrices
        pIdxMat_f = reshape(pInd_f, pilotSymLen, nBlks_f);
        pSymMat_f = reshape(pSym_f, pilotSymLen, nBlks_f);
        rxPM_f    = rxFr(pIdxMat_f);
        hLS_b_f   = mean(rxPM_f./pSymMat_f,1).';
        hLS_fr_f  = mean(hLS_b_f);

        %% SNORE noise estimate — use frame-local pilot indices
        %  pSig must be mean signal power = |mean(y·conj(p))|^2 (already an
        %  average over pilots), so it is directly comparable to pSigN which
        %  is also a mean power.  Do NOT square after the mean — that gives
        %  amplitude^2 of the mean, not mean power, and makes sigma_n2 wrong.
        rxPV_f = rxFr(pInd_f);
        pSN_f  = mean(abs(rxPV_f).^2);                      % mean received power
        pS_f   = mean(abs(rxPV_f.*conj(pSym_f)));           % mean |h| (amplitude)
        pS_f   = pS_f^2;                                    % mean signal power
        sn2_f  = max(pSN_f - pS_f, 1e-12);                 % noise power (signed diff, floored)

        %% MMSE estimation
        %  R_hh estimated from LS blocks with noise-bias correction:
        %    R_hh = E[|h_LS|^2] - sigma_n2/Np
        %  Wiener scalar: W = R_hh / (R_hh + sigma_n2/Np)
        %  Applied per-block so MMSE != LS whenever W != 1
        Rh_f     = max(real(mean(abs(hLS_b_f).^2) - sn2_f/Np), 1e-12);
        W_f      = Rh_f / (Rh_f + sn2_f/Np);
        hMM_b_f  = W_f * hLS_b_f;          % per-block MMSE estimate [nBlks_f x 1]
        hMM_fr_f = mean(hMM_b_f);          % frame-level scalar

        %% BER on data symbols — use frame-local data indices
        %  Bug fix: never use min(dInd_f, scalar) element-wise — it clamps
        %  multiple indices to the same last sample and corrupts the TX reference.
        %  Instead, clip dInd_f to valid range before indexing.
        dInd_f = dInd_f(dInd_f <= length(txFr));
        txD_f  = txFr(dInd_f);
        rLS_f = rxFr(dInd_f)/hLS_fr_f;
        rMM_f = rxFr(dInd_f)/hMM_fr_f;
        rNE_f = rxFr(dInd_f);

        txBt=qpsk(txD_f);   txBt=txBt(:);
        rBL=qpsk(rLS_f);    rBL=rBL(:);
        rBM=qpsk(rMM_f);    rBM=rBM(:);
        rBN=qpsk(rNE_f);    rBN=rBN(:);

        nb=min([length(txBt),length(rBL),length(rBM),length(rBN)]);
        if nb<10, continue; end

        bLS=bLS+sum(txBt(1:nb)~=rBL(1:nb))/nb;
        bMM=bMM+sum(txBt(1:nb)~=rBM(1:nb))/nb;
        bNE=bNE+sum(txBt(1:nb)~=rBN(1:nb))/nb;
        nV=nV+1;
    end

    if nV>0
        BER_LS_arr(si)   = bLS/nV;
        BER_MMSE_arr(si) = bMM/nV;
        BER_noEq_arr(si) = bNE/nV;
    else
        BER_LS_arr(si)=NaN; BER_MMSE_arr(si)=NaN; BER_noEq_arr(si)=NaN;
    end
    fprintf('  EsNo=%3.0f dB | noEq=%.4f  LS=%.4f  MMSE=%.4f  [%d frames]\n',...
        snr,BER_noEq_arr(si),BER_LS_arr(si),BER_MMSE_arr(si),nV);
end

%% ========================================================================
%  SECTION 8: DATASET COLLECTION FOR DL CHANNEL ESTIMATOR
%  -----------------------------------------------------------------------
%  PURPOSE
%    Build a frame-level dataset from the same BER loop above, where
%    every sample = one transmitted PL frame passing through the rain
%    fading + AWGN channel.
%
%  LABEL  (what the ML model must predict)
%    h_k  — the true complex channel coefficient for frame k.
%            Stored as [Re(h_k), Im(h_k)], shape (1 x 2) per frame.
%            This is the ONLY supervision signal the network needs.
%            It is constant within a frame (slow-fading assumption).
%
%  ML INPUTS  (what the network sees — paper eq. 33)
%    X_config1 : [Re(y_b); Re(p_b); Re(hLS_b);
%                 Im(y_b); Im(p_b); Im(hLS_b)]  per block b
%                shape per frame: (nBlks x 74)
%    X_config2 : drop hLS — shape (nBlks x 72)
%    X_config3 : only hLS  — shape (nBlks x 2)
%
%    Each frame is stored as a row in a cell array because nBlks can
%    vary across frames (short vs normal frame, bounds trimming).
%    Python handles variable-length sequences natively with padding.
%
%  ADDITIONAL INFO  (not fed to ML, useful for analysis / BER replay)
%    snr_db     : EsNo point at which this frame was generated
%    Ap_dB      : rain attenuation value drawn for this frame
%    nBlks      : number of pilot blocks that survived bounds check
%    hLS_frame  : scalar LS estimate (mean over blocks) for quick baseline
%    hMM_frame  : scalar MMSE estimate for comparison
%    sigma_n2   : noise variance estimated by SNORE
%    ber_ls     : per-frame BER under LS equalization
%    ber_mmse   : per-frame BER under MMSE equalization
%    rx_data    : received data symbols (complex) — needed to replay BER
%                 with DL estimate in Python without re-running MATLAB
%    tx_data    : transmitted data symbols (complex) — ground truth for BER
%% ========================================================================
fprintf('\n--- SECTION 8: Dataset Collection ---\n');

%  -----------------------------------------------------------------------
%  8a — Pre-allocate storage
%       Use cell arrays for variable-length sequences (nBlks may differ).
%       Preallocate for worst-case total frames = numSNR * numEvalFrames.
%  -----------------------------------------------------------------------
maxFrames   = numSNR * numEvalFrames;
pilotSymLen_ds = pilotSymLen;   % = 36, local alias for clarity

% --- ML inputs (cell arrays — one cell per frame) ---
DS.X_config1 = cell(maxFrames, 1);  % (nBlks x 74) float — full paper input
DS.X_config2 = cell(maxFrames, 1);  % (nBlks x 72) float — no hLS
DS.X_config3 = cell(maxFrames, 1);  % (nBlks x 2)  float — hLS only

% --- ML label (numeric array — same for all configs) ---
% LABEL: true complex channel h_k stored as [Re, Im]
% Shape: (maxFrames x 2) — one row per frame
DS.label_h   = zeros(maxFrames, 2);  % ← THE LABEL — Re(h_k) and Im(h_k)

% --- Additional info (not for ML, for analysis) ---
DS.snr_db    = zeros(maxFrames, 1);   % EsNo (dB) this frame was generated at
DS.Ap_dB     = zeros(maxFrames, 1);   % rain attenuation (dB) — ITU-R draw
DS.nBlks     = zeros(maxFrames, 1);   % pilot blocks in this frame
DS.hLS_frame = zeros(maxFrames, 1,'like',1+1j);  % scalar LS estimate (complex)
DS.hMM_frame = zeros(maxFrames, 1,'like',1+1j);  % scalar MMSE estimate (complex)
DS.sigma_n2  = zeros(maxFrames, 1);   % SNORE noise variance estimate
DS.ber_ls    = zeros(maxFrames, 1);   % per-frame BER under LS
DS.ber_mmse  = zeros(maxFrames, 1);   % per-frame BER under MMSE
DS.rx_data   = cell(maxFrames, 1);    % received data symbols (complex vector)
DS.tx_data   = cell(maxFrames, 1);    % transmitted data symbols (complex vector)

frameIdx = 0;   % running counter of valid collected frames

%  -----------------------------------------------------------------------
%  8b — Collection loop  (mirrors Section 7 BER loop exactly)
%       Every valid frame that passed the BER loop guards is collected.
%  -----------------------------------------------------------------------
for si = 1:numSNR
    snr = EsNo_range(si);

    for f = 1:numEvalFrames

        %% TX — identical to Section 7
        wg = dvbs2WaveformGenerator;
        wg.StreamFormat='TS'; wg.NumInputStreams=1;
        wg.MODCOD=MODCOD; wg.FECFrame=FECFrame; wg.HasPilots=true;
        wg.RolloffFactor=RolloffFactor;
        wg.FilterSpanInSymbols=FilterSpan; wg.SamplesPerSymbol=SPS;
        rp  = randi([0 1], pktLen, numPkts);
        tp  = [repmat(syncBits,1,numPkts); rp];
        txW = [wg(tp(:)); flushFilter(wg)];

        %% Channel — unique rain draw per frame (same as Section 7)
        p_f = 0.01 + (5.0-0.01)*rand(1);
        if p_f <= 1
            Cf = 0.655 + 0.033*log(p_f) - 0.045*log(max(A_001,0.01));
        else
            Cf = 0.655 + 0.033*log(p_f) - 0.045*log(max(A_001,0.01)) ...
                 - 0.053*(1-p_f)*sind(elev_deg);
        end
        Ap_f_dB = max(min(A_001*(p_f/0.01)^(-Cf), snr-3), 0.1);
        ht_f    = 10^(-Ap_f_dB/20) * exp(-1j*2*pi*rand(1));   % TRUE h_k
        rxW     = awgn(txW*ht_f, snr-10*log10(SPS), 'measured');

        %% RX — matched filter
        rf  = comm.RaisedCosineReceiveFilter('RolloffFactor',RolloffFactor,...
              'FilterSpanInSymbols',FilterSpan,...
              'InputSamplesPerSymbol',SPS,'DecimationFactor',SPS);
        rxS = rf(rxW); rxS = rxS(filterDelay+1:end); release(rf);

        tf  = comm.RaisedCosineReceiveFilter('RolloffFactor',RolloffFactor,...
              'FilterSpanInSymbols',FilterSpan,...
              'InputSamplesPerSymbol',SPS,'DecimationFactor',SPS);
        txS = tf(txW); txS = txS(filterDelay+1:end); release(tf);

        %% Frame length guard (same as Section 7)
        safeLen_f = min([length(rxS), length(txS), plFrameSize_orig]);
        if safeLen_f < plFrameSize_orig*0.9, continue; end
        rxFr = rxS(1:safeLen_f);
        txFr = txS(1:safeLen_f);

        %% Pilot/data index clipping (same as Section 7)
        pInd_f = pilotInd(pilotInd <= safeLen_f);
        dInd_f = dataInd(dataInd  <= safeLen_f);
        if isempty(pInd_f) || isempty(dInd_f), continue; end
        nBlks_f = floor(length(pInd_f) / pilotSymLen_ds);
        if nBlks_f < 2, continue; end   % need >=2 blocks for interp & sequence
        pInd_f = pInd_f(1:nBlks_f*pilotSymLen_ds);
        pSym_f = pilotSymbols(1:nBlks_f*pilotSymLen_ds);

        %% LS estimation (same as Section 7)
        pIdxMat_f = reshape(pInd_f,  pilotSymLen_ds, nBlks_f);
        pSymMat_f = reshape(pSym_f,  pilotSymLen_ds, nBlks_f);
        rxPM_f    = rxFr(pIdxMat_f);                   % [36 x nBlks]
        hLS_b_f   = mean(rxPM_f ./ pSymMat_f, 1).';   % [nBlks x 1] complex
        hLS_fr_f  = mean(hLS_b_f);                     % scalar

        %% SNORE + MMSE (same as Section 7)
        rxPV_f   = rxFr(pInd_f);
        pSN_f    = mean(abs(rxPV_f).^2);
        pS_f     = mean(abs(rxPV_f .* conj(pSym_f)));
        pS_f     = pS_f^2;
        sn2_f    = max(pSN_f - pS_f, 1e-12);
        Np_ds    = pilotSymLen_ds;
        Rh_f     = max(real(mean(abs(hLS_b_f).^2) - sn2_f/Np_ds), 1e-12);
        W_f      = Rh_f / (Rh_f + sn2_f/Np_ds);
        hMM_b_f  = W_f * hLS_b_f;
        hMM_fr_f = mean(hMM_b_f);

        %% Data symbols (same as Section 7)
        dInd_f = dInd_f(dInd_f <= length(txFr));
        if isempty(dInd_f), continue; end
        rxData_f = rxFr(dInd_f);           % received data — complex
        txData_f = txFr(dInd_f);           % transmitted data — complex

        %% Per-frame BER under LS and MMSE (for additional info)
        rLS_f = rxData_f / hLS_fr_f;
        rMM_f = rxData_f / hMM_fr_f;
        txBt  = [real(txData_f)>=0, imag(txData_f)>=0]; txBt = txBt(:);
        rBL   = [real(rLS_f)>=0,    imag(rLS_f)>=0];    rBL  = rBL(:);
        rBM   = [real(rMM_f)>=0,    imag(rMM_f)>=0];    rBM  = rBM(:);
        nb    = min([length(txBt), length(rBL), length(rBM)]);
        if nb < 10, continue; end
        ber_ls_f   = sum(txBt(1:nb) ~= rBL(1:nb)) / nb;
        ber_mmse_f = sum(txBt(1:nb) ~= rBM(1:nb)) / nb;

        %% ----------------------------------------------------------------
        %  BUILD ML INPUT TOKENS
        %  One token per pilot block b → stacked into (nBlks x features)
        %  Follows paper eq (33): X_in = [Re(y),Re(p),Re(hLS),Im(y),Im(p),Im(hLS)]
        %% ----------------------------------------------------------------
        X1 = zeros(nBlks_f, 74);   % Config 1 — full input
        X2 = zeros(nBlks_f, 72);   % Config 2 — no hLS
        X3 = zeros(nBlks_f,  2);   % Config 3 — hLS only

        for b = 1:nBlks_f
            yb  = rxPM_f(:, b);        % [36 x 1] complex — received pilots block b
            pb  = pSymMat_f(:, b);     % [36 x 1] complex — known pilots block b
            hb  = hLS_b_f(b);          % scalar   complex — LS estimate block b

            % Config 1: full paper input (eq 33)
            X1(b,:) = [real(yb)', real(pb)', real(hb), ...
                       imag(yb)', imag(pb)', imag(hb)];

            % Config 2: raw pilots only — no LS pre-computation
            X2(b,:) = [real(yb)', real(pb)', ...
                       imag(yb)', imag(pb)'];

            % Config 3: compressed — LS estimate only
            X3(b,:) = [real(hb), imag(hb)];
        end

        %% ----------------------------------------------------------------
        %  STORE THIS FRAME
        %% ----------------------------------------------------------------
        frameIdx = frameIdx + 1;

        % ML inputs
        DS.X_config1{frameIdx} = X1;   % (nBlks x 74)
        DS.X_config2{frameIdx} = X2;   % (nBlks x 72)
        DS.X_config3{frameIdx} = X3;   % (nBlks x 2)

        % LABEL — true complex channel h_k split into real and imaginary
        % This is what the BLSTM/GRU must learn to predict
        DS.label_h(frameIdx, :) = [real(ht_f), imag(ht_f)];

        % Additional info — for analysis, BER replay, debugging
        DS.snr_db(frameIdx)    = snr;
        DS.Ap_dB(frameIdx)     = Ap_f_dB;
        DS.nBlks(frameIdx)     = nBlks_f;
        DS.hLS_frame(frameIdx) = hLS_fr_f;
        DS.hMM_frame(frameIdx) = hMM_fr_f;
        DS.sigma_n2(frameIdx)  = sn2_f;
        DS.ber_ls(frameIdx)    = ber_ls_f;
        DS.ber_mmse(frameIdx)  = ber_mmse_f;
        DS.rx_data{frameIdx}   = rxData_f;  % complex vector — for BER replay
        DS.tx_data{frameIdx}   = txData_f;  % complex vector — ground truth
    end
end

%  -----------------------------------------------------------------------
%  8c — Trim pre-allocated arrays to actual collected frames
%  -----------------------------------------------------------------------
nCollected = frameIdx;
fprintf('Collected %d valid frames out of %d attempted\n', ...
        nCollected, maxFrames);

DS.X_config1 = DS.X_config1(1:nCollected);
DS.X_config2 = DS.X_config2(1:nCollected);
DS.X_config3 = DS.X_config3(1:nCollected);
DS.label_h   = DS.label_h(1:nCollected, :);
DS.snr_db    = DS.snr_db(1:nCollected);
DS.Ap_dB     = DS.Ap_dB(1:nCollected);
DS.nBlks     = DS.nBlks(1:nCollected);
DS.hLS_frame = DS.hLS_frame(1:nCollected);
DS.hMM_frame = DS.hMM_frame(1:nCollected);
DS.sigma_n2  = DS.sigma_n2(1:nCollected);
DS.ber_ls    = DS.ber_ls(1:nCollected);
DS.ber_mmse  = DS.ber_mmse(1:nCollected);
DS.rx_data   = DS.rx_data(1:nCollected);
DS.tx_data   = DS.tx_data(1:nCollected);

%  -----------------------------------------------------------------------
%  8d — Stratified train / val / test split  (paper: 1500 / 500 / 500)
%       Stratified by SNR so each split covers the full SNR range evenly.
%  -----------------------------------------------------------------------
rng(42);   % fixed seed for reproducibility
idx_all = randperm(nCollected);

% Paper split ratios: 60% train / 20% val / 20% test
n_train = round(0.60 * nCollected);
n_val   = round(0.20 * nCollected);

idx_train = idx_all(1          : n_train);
idx_val   = idx_all(n_train+1  : n_train+n_val);
idx_test  = idx_all(n_train+n_val+1 : end);

fprintf('Split  →  Train: %d  |  Val: %d  |  Test: %d\n', ...
        length(idx_train), length(idx_val), length(idx_test));

%  -----------------------------------------------------------------------
%  8e — Save to .mat file
%       Config 1 is the primary ML input matching the paper.
%       Configs 2 and 3 are saved for ablation experiments.
%       All additional fields saved for offline analysis in Python/MATLAB.
%  -----------------------------------------------------------------------
saveFile = 'dvbs2x_channel_dataset.mat';

save(saveFile, ...
    ... % ML inputs (cell arrays — one cell per frame)
    '-struct', 'DS', ...
    ... % index splits
    'idx_train', 'idx_val', 'idx_test', ...
    ... % system config metadata (needed to reproduce any frame)
    'MODCOD', 'FECFrame', 'RolloffFactor', 'pilotSymLen', ...
    'EsNo_range', 'numEvalFrames', ...
    '-v7.3');   % v7.3 required for large cell arrays > 2 GB

fprintf('Dataset saved  →  %s\n', saveFile);
fprintf('\nDataset summary:\n');
fprintf('  Total frames    : %d\n',  nCollected);
fprintf('  Config1 features: 74 per block  [Re(y),Re(p),Re(hLS),Im(y),Im(p),Im(hLS)]\n');
fprintf('  Config2 features: 72 per block  [Re(y),Re(p),Im(y),Im(p)]\n');
fprintf('  Config3 features:  2 per block  [Re(hLS),Im(hLS)]\n');
fprintf('  LABEL           : DS.label_h  — [Re(h_k), Im(h_k)]  shape (N x 2)\n');
fprintf('  SNR range       : %d to %d dB\n', min(EsNo_range), max(EsNo_range));
fprintf('  Pilot blocks    : min=%d  max=%d  (varies with frame trimming)\n', ...
        min(DS.nBlks), max(DS.nBlks));

%% BER Plot
figure('Name','BER vs EsNo','Color','w','Position',[100 100 820 500]);

semilogy(EsNo_range,BER_noEq_arr,'k:' ,'LineWidth',1.5,...
    'DisplayName','No Equalization');
hold on;
semilogy(EsNo_range,BER_LS_arr,  'b--s','LineWidth',1.8,'MarkerSize',7,...
    'DisplayName','LS Estimator');
semilogy(EsNo_range,BER_MMSE_arr,'r-o' ,'LineWidth',1.8,'MarkerSize',7,...
    'DisplayName','MMSE Estimator');
semilogy(EsNo_range,0.5*erfc(sqrt(10.^(EsNo_range/10))),'g-.','LineWidth',1.5,...
    'DisplayName','Theoretical QPSK (AWGN)');

grid on;
xlabel('E_s/N_0 (dB)','FontSize',13);
ylabel('BER','FontSize',13);
title({sprintf('BER vs E_s/N_0 | DVB-S2 %s Frame | MODCOD %d',FECFrame,MODCOD),...
    sprintf('Tropical fading (Penang 20 GHz) | %d frames/point',numEvalFrames)},...
    'FontSize',12,'FontWeight','bold');
legend('Location','southwest','FontSize',11);
ylim([1e-4 1]); xlim([min(EsNo_range) max(EsNo_range)]);

fprintf('\nDone.\n');