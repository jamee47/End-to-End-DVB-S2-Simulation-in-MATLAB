clc; clearvars; close all;

% ══════════════════════════════════════════════════════════════════════════
%  DVB-S2 Synthetic Dataset Generator
%  Channel Model : ITU-R P.618 Rain Attenuation (from experimental.m)
%  Output Schema : real_pilots | imag_pilots | htrue_real | htrue_imag |
%                  snr_dB | nvar | modcod | rainAtt_dB
%
%  Effective Es/N0 at receiver = esno_dB - rainAtt_dB
%  (rain attenuation is subtracted from the nominal link Es/N0)
%
%  Notes
%  ─────
%  • Pilots are stored as a flat row: [P×1 real part, P×1 imag part]
%    where P = numPilotBlks × pilotSymLen (all pilot symbols in the frame).
%  • H_true is the full composite channel coefficient (rain + free-space +
%    beam gain) exactly as computed in experimental.m (Eq. 7–8 region).
%  • snr_dB stored here is the *effective* post-rain Es/N0 in dB:
%        snr_dB = esno_dB - rainAtt_dB
%  • nvar is estimated from received pilots after matched filtering.
%
%  Valid MODCOD IDs (normal FEC frame, dvbs2WaveformGenerator compatible)
%  ───────────────────────────────────────────────────────────────────────
%   ID │ Modulation │ Code Rate │ Min Es/No (dB)
%  ────┼────────────┼───────────┼───────────────
%    1 │ QPSK       │ 1/4       │  -2.4
%    2 │ QPSK       │ 1/2       │   0.0
%    3 │ QPSK       │ 3/4       │   1.8
%    4 │ 8PSK       │ 2/3       │   4.0
%    5 │ 8PSK       │ 3/4       │   5.0
%    6 │ 8PSK       │ 5/6       │   6.0
%    7 │ 16APSK     │ 2/3       │   7.5
%    8 │ 16APSK     │ 3/4       │   8.8
%    9 │ 16APSK     │ 5/6       │   9.8
%   10 │ 32APSK     │ 3/4       │  12.0
%   11 │ 32APSK     │ 5/6       │  13.5
%   12 │ 32APSK     │ 8/9       │  15.0
%   13 │ 32APSK     │ 9/10      │  16.0
%  ────┴────────────┴───────────┴───────────────
%  Single MODCOD  →  modcod_selection = 5;
%  Multiple       →  modcod_selection = [2 5 8 11];
%  All valid IDs  →  modcod_selection = 1:13;
% ══════════════════════════════════════════════════════════════════════════

%% ═══════════════════════════════════════════════════════════════════════
%                        USER CONFIGURATION
%  ▼ Edit the values in this section only ▼
%% ═══════════════════════════════════════════════════════════════════════

% ── MODCOD selection ──────────────────────────────────────────────────
%  Specify one or more MODCOD IDs from the table above.
%  The dataset will cycle randomly across all IDs listed here.
modcod_selection = [2 3 5 6 11 12];

% ── Number of frames and output ───────────────────────────────────────
numFrames = 50;            % total frames to generate
saveEvery = 5;             % MAT checkpoint interval (frames)
outDir    = '../dataset_output';

% ── SNR diversity (nominal Es/N0 before rain, in dB) ─────────────────
esno_pool = 15;

% ── Temporal rain-fade model ──────────────────────────────────────────
%  Rain attenuation evolves as a 3-state Markov chain:
%    CLEAR  : rainAtt_dB ~ 0–1 dB   (clear sky / light cloud)
%    ONSET  : rainAtt_dB ~ 1–6 dB   (rain building up)
%    FADE   : rainAtt_dB ~ 6–15 dB  (deep fade event)
%
%  Transition probabilities per frame (tune to change fade duration):
%    stay_clear  — probability of remaining in CLEAR  (high → long clear runs)
%    stay_onset  — probability of remaining in ONSET  (controls ramp duration)
%    stay_fade   — probability of remaining in FADE   (high → sustained fades)
%
%  Ka-band 20 GHz tropical reference: deep fades up to ~15 dB are realistic
%  for a worst-case 0.01% exceedance. Values beyond 15 dB are rare outages.

fade_stay_clear = 0.92;   % prob of staying CLEAR  each frame
fade_stay_onset = 0.60;   % prob of staying ONSET  each frame
fade_stay_fade  = 0.80;   % prob of staying FADE   each frame

%% ═══════════════════════════════════════════════════════════════════════
%                     END OF USER CONFIGURATION
%% ═══════════════════════════════════════════════════════════════════════

%% ── Validate MODCOD selection ────────────────────────────────────────────
valid_ids = 1:13;
assert(all(ismember(modcod_selection, valid_ids)), ...
    'Invalid MODCOD ID detected. Valid range is 1–13 (normal FEC frame).');
modcod_pool = modcod_selection(:)';   % row vector, alias used throughout

FECframe = 'normal';          % FEC frame type — must be set before filename tag

%% ── Output files (filename encodes modulation order + code rate) ─────────
%  Each selected MODCOD contributes one token: <Modulation>_<num>_<den>
%  e.g. MODCOD 5 (8PSK 3/4)  → "8PSK_3_4"
%       MODCOD 2 (QPSK 1/2)  → "QPSK_1_2"
%  Multiple MODCODs joined with double underscores:
%  e.g. [2 5 11] → "QPSK_1_2__8PSK_3_4__32APSK_5_6"
modName = containers.Map([4 8 16 32], {'QPSK','8PSK','16APSK','32APSK'});

mc_tags = cell(1, numel(modcod_pool));
for ii = 1 : numel(modcod_pool)
    [mOrd, cRate] = satcom.internal.dvbs.getS2PHYParams(modcod_pool(ii), FECframe);
    [cr_num, cr_den] = rat(cRate, 1e-4);
    mc_tags{ii} = sprintf('%s_%d_%d', modName(mOrd), cr_num, cr_den);
end
mc_tag = strjoin(mc_tags, '__');

if ~exist(outDir, 'dir'), mkdir(outDir); end
matFile = fullfile(outDir, sprintf('dvbs2_dataset_%s.mat', mc_tag));
csvFile = fullfile(outDir, sprintf('dvbs2_dataset_%s.csv', mc_tag));

%% ── Tropical channel parameters (from experimental.m) ────────────────────
chanParams.Frequency_GHz    = 20;
chanParams.Latitude         = 5.17;
chanParams.Altitude_m       = 57;
chanParams.ElevationAngle   = 45;
chanParams.Polarization     = 'V';
chanParams.R001_mmph        = 130;
chanParams.IsothermHeight_km = 4.5;

%% ── Free-space link constants ────────────────────────────────────────────
c_light   = 3e8;
kb        = 1.38064852e-23;
freq_hz   = chanParams.Frequency_GHz * 1e9;
lambda    = c_light / freq_hz;
d0        = 35788e3;

simParams.G_s_dBi      = 52.0;
simParams.G_r_dBi      = 41.7;
simParams.NoiseBW_Hz   = 50e6;
simParams.NoiseTemp_K  = 207;
simParams.Theta3dB_deg = 0.4;
simParams.BeamOffset_deg = 0.0;

Gr_linear = 10^(simParams.G_r_dBi / 10);
b_max     = (lambda / (4*pi))^2 * (1/(d0^2)) * ...
            (Gr_linear / (kb * simParams.NoiseBW_Hz * simParams.NoiseTemp_K));
Gs_linear = 10^(simParams.G_s_dBi / 10);
theta3dB  = deg2rad(simParams.Theta3dB_deg);
theta_off = deg2rad(simParams.BeamOffset_deg);
u_val     = 2.07123 * sin(theta_off) / sin(theta3dB);
if abs(u_val) < 1e-8
    beamPattern = 1;
else
    beamPattern = besselj(1,u_val)/(2*u_val) + 36*besselj(3,u_val)/(u_val^3);
end
b_gain = Gs_linear * abs(beamPattern)^2;

%% ── Waveform base config ─────────────────────────────────────────────────
RolloffFactor       = 0.35;
FilterSpanInSymbols = 10;
SPS                 = 2;
HasPilots           = true;
filterDelay         = FilterSpanInSymbols / 2;
syncBits            = [0;1;0;0;0;1;1;1];
pktLen              = 1496;

%% ── Pre-compute PHY geometry per MODCOD ──────────────────────────────────
%  Use containers.Map keyed by MODCOD id to avoid struct-array index gaps
%  when the selected MODCODs are non-contiguous (e.g. [2 5 11]).
geom = containers.Map('KeyType','int32','ValueType','any');

for mc = modcod_pool
    [modOrd, cRate, cwLen] = satcom.internal.dvbs.getS2PHYParams(mc, FECframe);
    dLen        = cwLen / log2(modOrd);
    slotLen     = 90;
    pbFreq      = 16;
    pilotSymLen = 36;
    nPBlks      = floor(dLen / (slotLen * pbFreq));
    if floor(dLen/(slotLen*pbFreq)) == dLen/(slotLen*pbFreq)
        nPBlks = nPBlks - 1;
    end
    pilotLen  = nPBlks * pilotSymLen;
    frmSz     = dLen + pilotLen + slotLen;

    plScr  = satcom.internal.dvbs.plScramblingIntegerSequence(0);
    cMap   = [1; 1j; -1; -1j];
    cSeq   = cMap(plScr + 1);

    [~, pInd_raw] = satcom.internal.dvbs.pilotBlock(nPBlks);
    pInd_raw = pInd_raw(:)';
    pInd     = pInd_raw + slotLen;
    refP     = (1+1j)/sqrt(2) .* cSeq(pInd_raw(:));

    g.modOrder       = modOrd;
    g.codeRate       = cRate;
    g.plFrameSize    = frmSz;
    g.numPilotBlks   = nPBlks;
    g.pilotSymLen    = pilotSymLen;
    g.pilotInd       = pInd;
    g.pilotInd_raw   = pInd_raw;
    g.refPilots      = refP;
    g.Np             = length(pInd);   % total pilot symbols per frame
    g.numPktsPerFrame = [];            % filled after wavegen init below
    geom(int32(mc)) = g;
end

%% ── Initialise wavegen to fetch MinNumPackets per MODCOD ─────────────────
wavegen = dvbs2WaveformGenerator;
wavegen.StreamFormat        = 'TS';
wavegen.NumInputStreams      = 1;
wavegen.FECFrame            = FECframe;
wavegen.HasPilots           = HasPilots;
wavegen.RolloffFactor       = RolloffFactor;
wavegen.FilterSpanInSymbols = FilterSpanInSymbols;
wavegen.SamplesPerSymbol    = SPS;

for mc = modcod_pool
    wavegen.MODCOD = mc;
    wavegen.DFL    = getDFL(mc, FECframe);
    g = geom(int32(mc));
    g.numPktsPerFrame = wavegen.MinNumPackets(1);
    geom(int32(mc)) = g;
end

%% ── Determine max Np across MODCODs for fixed-width CSV ──────────────────
maxNp = max(cellfun(@(mc) geom(int32(mc)).Np, num2cell(modcod_pool)));

%% ── Open CSV and write header ────────────────────────────────────────────
fid = fopen(csvFile, 'w');
if fid == -1, error('Cannot open CSV: %s', csvFile); end

hdr = 'frame_idx,modcod,snr_dB,nvar,rainAtt_dB,htrue_real,htrue_imag,';
for p = 1:maxNp, hdr = [hdr, sprintf('real_pilot_%d,',  p)]; end %#ok<AGROW>
for p = 1:maxNp, hdr = [hdr, sprintf('imag_pilot_%d,',  p)]; end %#ok<AGROW>
hdr(end) = [];   % remove trailing comma
fprintf(fid, '%s\n', hdr);
fprintf('CSV header written (%d pilot columns).\n', maxNp*2);

%% ── Pre-allocate dataset struct ──────────────────────────────────────────
dataset(numFrames) = struct( ...
    'frame_idx',   [], ...   % scalar
    'modcod',      [], ...   % MODCOD index (integer)
    'snr_dB',      [], ...   % effective Es/N0 = esno_dB - rainAtt_dB  [dB]
    'nvar',        [], ...   % noise variance estimated from pilots
    'rainAtt_dB',  [], ...   % ITU-R P.618 rain attenuation             [dB]
    'htrue_real',  [], ...   % real(H_true) — composite channel coeff
    'htrue_imag',  [], ...   % imag(H_true)
    'real_pilots', [], ...   % [Np × 1] real part of received pilot symbols
    'imag_pilots', []);      % [Np × 1] imag part of received pilot symbols

%% ── Progress header ──────────────────────────────────────────────────────
fprintf('\nGenerating %d frames …\n\n', numFrames);
fprintf('%-7s %-8s %-12s %-12s %-10s\n', ...
    'Frame', 'MODCOD', 'Es/No(nom)', 'RainAtt(dB)', 'SNR_eff(dB)');
fprintf('%s\n', repmat('-', 1, 55));

%% ── Initialise temporal fade state ──────────────────────────────────────
fadeState = 1;   % start in CLEAR (1=CLEAR, 2=ONSET, 3=FADE)

%% ── Main generation loop ─────────────────────────────────────────────────
for frameIdx = 1 : numFrames

    % ── Random selections ─────────────────────────────────────────────────
    mc       = modcod_pool(randi(length(modcod_pool)));
    esno_dB  = esno_pool(randi(length(esno_pool)));
    g        = geom(int32(mc));

    % ── Temporal rain-fade state machine ──────────────────────────────────
    %  States: 1 = CLEAR, 2 = ONSET, 3 = FADE
    %  Transitions evaluated once per frame; rainAtt_dB sampled from the
    %  uniform range of the current state.
    u = rand();
    switch fadeState
        case 1  % CLEAR
            if u < fade_stay_clear
                fadeState = 1;
            else
                fadeState = 2;          % move to onset
            end
        case 2  % ONSET
            if u < fade_stay_onset
                fadeState = 2;
            elseif u < fade_stay_onset + (1-fade_stay_onset)*0.6
                fadeState = 3;          % deepen to full fade
            else
                fadeState = 1;          % dissipate back to clear
            end
        case 3  % FADE
            if u < fade_stay_fade
                fadeState = 3;
            else
                fadeState = 2;          % begin recovery through onset
            end
    end

    switch fadeState
        case 1,  rainAtt_dB = 0   + 1   * rand();   % 0–1  dB  clear sky
        case 2,  rainAtt_dB = 1   + 5   * rand();   % 1–6  dB  onset/recovery
        case 3,  rainAtt_dB = 6   + 9   * rand();   % 6–15 dB  deep fade
    end

    % ── Effective SNR = nominal Es/N0 minus rain loss ─────────────────────
    %    This is the true quantity of interest for the learning task.
    snr_eff_dB = esno_dB - rainAtt_dB;
    snr_awgn   = snr_eff_dB - 10*log10(SPS);   % matched-filter SNR

    % ── Build composite channel coefficient (mirrors experimental.m) ──────
    %    Ap = amplitude attenuation factor (> 1 means loss)
    Ap             = 10^(rainAtt_dB / 20);
    phi            = 2 * pi * rand(1);
    h_attenuation  = (1 / sqrt(Ap)) * exp(-1i * phi);
    H_true         = h_attenuation * sqrt(b_gain) * sqrt(b_max);

    % ── Transmit waveform ─────────────────────────────────────────────────
    release(wavegen);
    wavegen.MODCOD = mc;
    wavegen.DFL    = getDFL(mc, FECframe);
    numPkts        = g.numPktsPerFrame;

    pn = comm.PNSequence('Polynomial', 'x9+x5+1', ...
        'InitialConditions', [zeros(1,8) 1], ...
        'VariableSizeOutput', true, ...
        'MaximumOutputSize', [pktLen * numPkts, 1]);
    reset(pn);
    rawBits   = pn(pktLen * numPkts);
    txRawPkts = reshape(rawBits, pktLen, numPkts);
    txPkts    = [repmat(syncBits, 1, numPkts); txRawPkts];
    inputBits = txPkts(:);

    txWaveform = [wavegen({inputBits}); flushFilter(wavegen)];

    % ── Apply channel: rain attenuation + phase rotation + AWGN ──────────
    attWaveform = H_true * txWaveform;
    rxNoisy     = awgn(attWaveform, snr_awgn, 'measured');
    clearvars txWaveform attWaveform rawBits txRawPkts txPkts inputBits pn;

    % ── Matched filter + symbol timing ────────────────────────────────────
    rxFilter = comm.RaisedCosineReceiveFilter( ...
        'RolloffFactor',         RolloffFactor, ...
        'FilterSpanInSymbols',   FilterSpanInSymbols, ...
        'InputSamplesPerSymbol', SPS, ...
        'DecimationFactor',      SPS);
    rxSym_raw  = rxFilter(rxNoisy);
    release(rxFilter);
    clearvars rxNoisy;

    % Symmetric truncation (as in experimental.m "FINE-TUNED" comment)
    rxSymbols = rxSym_raw(filterDelay + 1 : end - filterDelay);
    clearvars rxSym_raw;

    if length(rxSymbols) < g.plFrameSize
        warning('Frame %d: insufficient symbols — skipping.', frameIdx);
        clearvars rxSymbols;
        continue;
    end

    % ── Extract first PL frame ────────────────────────────────────────────
    rxFrame    = rxSymbols(1 : g.plFrameSize);
    clearvars rxSymbols;

    % ── Received pilot symbols ────────────────────────────────────────────
    rxPilotVec = rxFrame(g.pilotInd(:));   % [Np × 1] complex

    % ── Noise variance estimate (pilot-based) ─────────────────────────────
    nvar = DVBS2NoiseVarEstimate(rxFrame, g.pilotInd, g.refPilots, false);
    clearvars rxFrame;

    % ── Store in struct ───────────────────────────────────────────────────
    dataset(frameIdx).frame_idx   = frameIdx;
    dataset(frameIdx).modcod      = mc;
    dataset(frameIdx).snr_dB      = snr_eff_dB;        % effective Es/N0
    dataset(frameIdx).nvar        = nvar;
    dataset(frameIdx).rainAtt_dB  = rainAtt_dB;
    dataset(frameIdx).htrue_real  = real(H_true);
    dataset(frameIdx).htrue_imag  = imag(H_true);
    dataset(frameIdx).real_pilots = real(rxPilotVec);   % [Np × 1]
    dataset(frameIdx).imag_pilots = imag(rxPilotVec);   % [Np × 1]

    % ── Write CSV row (pad pilots to maxNp with NaN) ──────────────────────
    rp_pad = NaN(1, maxNp);  ip_pad = NaN(1, maxNp);
    rp_pad(1:g.Np) = real(rxPilotVec)';
    ip_pad(1:g.Np) = imag(rxPilotVec)';

    rowVec = [frameIdx, mc, snr_eff_dB, nvar, rainAtt_dB, ...
              real(H_true), imag(H_true), rp_pad, ip_pad];
    fprintf(fid, '%s\n', strjoin(arrayfun(@(v) sprintf('%.8g', v), ...
        rowVec, 'UniformOutput', false), ','));

    % ── Console progress ──────────────────────────────────────────────────
    fprintf('%-7d %-8d %-12.2f %-12.4f %-10.2f\n', ...
        frameIdx, mc, esno_dB, rainAtt_dB, snr_eff_dB);

    % ── Cleanup ───────────────────────────────────────────────────────────
    clearvars rxPilotVec rxFrame nvar H_true h_attenuation Ap phi ...
        rp_pad ip_pad rowVec mc esno_dB rainAtt_dB snr_eff_dB ...
        snr_awgn g numPkts;

    % ── MAT checkpoint ────────────────────────────────────────────────────
    if mod(frameIdx, saveEvery) == 0 || frameIdx == numFrames
        save(matFile, 'dataset', '-v7.3');
        fprintf('  → Checkpoint %d/%d saved → %s\n', frameIdx, numFrames, matFile);
    end
end

fclose(fid);

fprintf('\n══ Dataset generation complete ══\n');
fprintf('  Frames : %d\n', numFrames);
fprintf('  MAT    → %s\n', matFile);
fprintf('  CSV    → %s\n', csvFile);

%% ══════════════════════════════════════════════════════════════════════════
%  Helper functions
%% ══════════════════════════════════════════════════════════════════════════

function dfl = getDFL(modCod, fecFrame)
%GETDFL Data field length for MODCOD + FEC frame type.
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
%DVBS2NOISEVARESTIMATE Pilot-based noise variance estimator.
if normFlag
    rxData = rxData / sqrt(mean(abs(rxData).^2));
end
rxPilots = rxData(pilotInd(:));
pSigN    = mean(abs(rxPilots).^2);
pSig     = abs(mean(rxPilots .* conj(refPilots(:)))).^2;
nVarEst  = abs(pSigN - pSig);
end


function Ap_dB = calcRainAttenuationP618(chanParams, p)
%CALCRAINSATTENUATIONP618  ITU-R P.618 rain attenuation (from experimental.m).
f_GHz  = chanParams.Frequency_GHz;
theta_deg = chanParams.ElevationAngle;
lat_deg   = chanParams.Latitude;
hs_km     = chanParams.Altitude_m / 1000;
R001      = chanParams.R001_mmph;
h0_km     = chanParams.IsothermHeight_km;
hR_km     = h0_km + 0.36;
Re_km     = 8500;

if R001 <= 0 || hR_km <= hs_km
    Ap_dB = 0; return;
end

if theta_deg >= 5
    Ls = (hR_km - hs_km) / sind(theta_deg);
else
    Ls = 2*(hR_km - hs_km) / (sqrt(sind(theta_deg)^2 + ...
         2*(hR_km - hs_km)/Re_km) + sind(theta_deg));
end
LG = Ls * cosd(theta_deg);

[k, alpha] = ituP838_k_alpha(f_GHz, theta_deg, chanParams.Polarization);
gamma_R    = k * (R001^alpha);

r001 = 1 / (1 + 0.78*sqrt((LG*gamma_R)/f_GHz) - 0.38*(1 - exp(-2*LG)));
zeta = atand((hR_km - hs_km) / (LG*r001));
if zeta > theta_deg
    LR = (LG*r001) / cosd(theta_deg);
else
    LR = (hR_km - hs_km) / sind(theta_deg);
end

chi = max(0, 36 - abs(lat_deg));
v001 = 1 / (1 + sqrt(sind(theta_deg)) * (31*(1 - exp(-theta_deg/(1+chi))) * ...
       (sqrt(LR*gamma_R)/(f_GHz^2)) - 0.45));
LE   = LR * v001;
A001 = gamma_R * LE;

if abs(p - 0.01) < 1e-12
    Ap_dB = A001; return;
end

if p >= 1 || abs(lat_deg) >= 36
    beta = 0;
elseif theta_deg >= 25
    beta = -0.005*(abs(lat_deg) - 36);
else
    beta = -0.005*(abs(lat_deg) - 36) + 1.8 - 4.25*sind(theta_deg);
end

if A001 <= 0
    Ap_dB = 0; return;
end

exponent = -(0.655 + 0.033*log(p) - 0.045*log(A001) - beta*(1-p)*sind(theta_deg));
Ap_dB    = A001 * (p/0.01)^exponent;
end


function [k, alpha] = ituP838_k_alpha(f_GHz, theta_deg, polarization)
%ITUP838_K_ALPHA  Specific attenuation coefficients (ITU-R P.838).
x = log10(f_GHz);

akH = [-5.33980, -0.35351, -0.23789, -0.94158];
bkH = [-0.10008,  1.26970,  0.86036,  0.64552];
ckH = [ 1.13098,  0.45400,  0.15354,  0.16817];
akV = [-3.80595, -3.44965, -0.39902,  0.50167];
bkV = [ 0.56934, -0.22911,  0.73042,  1.07319];
ckV = [ 0.81061,  0.51059,  0.11899,  0.27195];

logkH = sum(akH .* exp(-((x-bkH)./ckH).^2)) - 0.18961*x + 0.71147;
logkV = sum(akV .* exp(-((x-bkV)./ckV).^2)) - 0.16398*x + 0.63297;
kH = 10^logkH;  kV = 10^logkV;

aaH = [-0.14318,  0.29591,  0.32177, -5.37610, 16.1721];
baH = [ 1.82442,  0.77564,  0.63773, -0.96230, -3.29980];
caH = [-0.55187,  0.19822,  0.13164,  1.47828,  3.43990];
aaV = [-0.07771,  0.56727, -0.20238, -48.2991, 48.5833];
baV = [ 2.33840,  0.95545,  1.14520,  0.791669, 0.791459];
caV = [-0.76284,  0.54039,  0.26809,  0.116226, 0.116479];

alphaH = sum(aaH .* exp(-((x-baH)./caH).^2)) + 0.67849*x - 1.95537;
alphaV = sum(aaV .* exp(-((x-baV)./caV).^2)) - 0.053739*x + 0.83433;

pol = upper(string(polarization));
if pol == "H",     tau_deg = 0;
elseif pol == "V", tau_deg = 90;
else,              tau_deg = 45;
end

k     = (kH + kV + (kH-kV)*cosd(theta_deg)^2*cosd(2*tau_deg)) / 2;
alpha = (kH*alphaH + kV*alphaV + (kH*alphaH - kV*alphaV)*cosd(theta_deg)^2*cosd(2*tau_deg)) / (2*k);
end