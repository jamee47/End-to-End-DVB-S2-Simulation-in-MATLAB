clc; clearvars; close all;

% ══════════════════════════════════════════════════════════════════════════════
%  DVB-S2 Dataset Generator  —  Paper-aligned version
%  Awad et al., "End-to-end DVB-S2X system design with DL-based channel
%  estimation over satellite fading channels at Ka-band", Comp. Networks 2023.
%
%  STRUCTURE (each section is self-contained and clearly labelled):
%  ─────────────────────────────────────────────────────────────────
%  SEC 1  │ USER CONFIGURATION          ← edit only here
%  SEC 2  │ ITU-R P.618 CHANNEL         ← tropical rain attenuation
%  SEC 3  │ LINK BUDGET & CHANNEL COEFF ← free-space, beam, h_k build
%  SEC 4  │ DVB-S2 PHY GEOMETRY         ← frame/pilot structure per MODCOD
%  SEC 5  │ OUTPUT FILE SETUP           ← MAT + CSV initialisation
%  SEC 6  │ TRANSMITTER                 ← waveform generation per frame
%  SEC 7  │ CHANNEL                     ← apply h_k + AWGN
%  SEC 8  │ RECEIVER                    ← matched filter, symbol timing
%  SEC 9  │ Xin ASSEMBLY & SAVING       ← paper Xin (no LS), MAT, CSV
%  ─────────────────────────────────────────────────────────────────
%
%  Xin saved (paper Eq. 33, LS part omitted as requested):
%      Xin = [ real(y_full)^T,  real(p_ref)^T,
%              imag(y_full)^T,  imag(p_ref)^T ]^T
%  where y_full = full received PLFRAME symbols (data + pilots, post-filter)
%        p_ref  = known reference pilot symbols (zero-padded to frame length)
%
%  CSV schema (one row per frame):
%      frame_idx | modcod | esno_nom_dB | rainAtt_dB | snr_eff_dB |
%      nvar | htrue_real | htrue_imag | Ap_linear | phi_rad |
%      b_gain | b_max | channel_type |
%      Xin_1 … Xin_N   (flattened, N = 4 × plFrameSize)
%
%  NOTE:  Full-frame CSV rows are very wide (4 × ~32 400 columns for QPSK
%         normal frame).  For large numFrames consider MAT-only storage and
%         set SAVE_CSV = false in SEC 1.
% ══════════════════════════════════════════════════════════════════════════════


%% ═══════════════════════════════════════════════════════════════════════════
%  SEC 1 — USER CONFIGURATION
%  ▼ Edit the values in this section only ▼
%% ═══════════════════════════════════════════════════════════════════════════

% ── MODCOD selection ─────────────────────────────────────────────────────────
%  Pick from the table below.  The loop randomly cycles across all listed IDs.
%
%   ID │ Modulation │ Code Rate │ Approx min Es/No
%  ────┼────────────┼───────────┼─────────────────
%    1 │ QPSK       │ 1/4       │  -2.4 dB
%    2 │ QPSK       │ 1/2       │   0.0 dB   ← paper QPSK proxy
%    3 │ QPSK       │ 3/4       │   1.8 dB
%    4 │ 8PSK       │ 2/3       │   4.0 dB
%    5 │ 8PSK       │ 3/4       │   5.0 dB
%    6 │ 8PSK       │ 5/6       │   6.0 dB
%    7 │ 16APSK     │ 2/3       │   7.5 dB
%    8 │ 16APSK     │ 3/4       │   8.8 dB
%    9 │ 16APSK     │ 5/6       │   9.8 dB
%   10 │ 32APSK     │ 3/4       │  12.0 dB   ← paper 32APSK proxy
%   11 │ 32APSK     │ 5/6       │  13.5 dB
%   12 │ 32APSK     │ 8/9       │  15.0 dB
%   13 │ 32APSK     │ 9/10      │  16.0 dB
modcod_selection = [2 5 10];      % IDs to include in dataset

% ── Frame count & checkpointing ──────────────────────────────────────────────
numFrames   = 1;       % total frames to generate (paper: 1500 train + 500 val + 500 test)
saveEvery   = 1;        % MAT checkpoint every N frames
SAVE_CSV    = true;       % set false to skip CSV (much faster for large runs)
outDir      = '../dataset_output';

% ── Nominal Es/No pool (dB) — BEFORE rain attenuation ────────────────────────
%  Each frame randomly draws one value.  Rain is then subtracted to get the
%  effective Es/No the receiver actually sees.
%  Rationale: sweeping this range ensures the DL model trains across the full
%  SNR range shown in the paper's BER/NMSE figures (≈5–35 dB effective).
esno_pool = 10 : 5 : 35;    % [10 15 20 25 30 35] dB  — edit freely

% ── FEC frame type ───────────────────────────────────────────────────────────
FECframe = 'normal';      % 'normal' (64 800 bits) | 'short' (16 200 bits)

% ── Waveform filter settings ─────────────────────────────────────────────────
RolloffFactor       = 0.20;   % paper uses 20 % roll-off (α = 0.20)
FilterSpanInSymbols = 10;
SPS                 = 2;      % samples per symbol

% ── Random seed (set to fixed integer for reproducibility, 0 = random) ───────
RNG_SEED = 0;


%% ═══════════════════════════════════════════════════════════════════════════
%  END OF USER CONFIGURATION — do not edit below unless you know what you do
%% ═══════════════════════════════════════════════════════════════════════════

if RNG_SEED ~= 0, rng(RNG_SEED); end


%% ═══════════════════════════════════════════════════════════════════════════
%  SEC 2 — ITU-R P.618 TROPICAL CHANNEL PARAMETERS
%  Reference: Penang, Malaysia (Table 1 of the paper)
%  ITU-R recommendations used: P.618, P.839, P.838
% ═══════════════════════════════════════════════════════════════════════════

% ── Site parameters (paper Table 1, tropical channel) ────────────────────────
itu.channel_type        = 'tropical';
itu.site                = 'Penang, Malaysia';
itu.Frequency_GHz       = 20;          % Ka-band downlink
itu.Latitude_deg        = 5.17;        % deg N
itu.Longitude_deg       = 100.4;       % deg E
itu.Altitude_m          = 57;          % m above sea level
itu.ElevationAngle_deg  = 45;          % satellite elevation angle
itu.Polarization        = 'V';         % vertical polarisation
itu.R001_mmph           = 130;         % rainfall rate exceeded 0.01 % of year
itu.IsothermHeight_km   = 4.5;         % 0°C isotherm height

% ── Compute ITU-R P.618 rain attenuation at p = 0.01 % (worst-case) ──────────
%  This gives A001 — the single deterministic reference attenuation used to
%  anchor the statistical distribution drawn each frame.
A001_dB = calcRainAttenuationP618(itu, 0.01);

fprintf('═══════════════════════════════════════════════════════\n');
fprintf('  ITU-R P.618 Rain Attenuation (tropical, 0.01%%)\n');
fprintf('  Site          : %s\n',   itu.site);
fprintf('  Frequency     : %.0f GHz\n', itu.Frequency_GHz);
fprintf('  Elevation     : %.0f deg\n', itu.ElevationAngle_deg);
fprintf('  R_0.01        : %.0f mm/h\n', itu.R001_mmph);
fprintf('  A_001 (ref)   : %.2f dB   ← worst-case 0.01%% exceedance\n', A001_dB);
fprintf('═══════════════════════════════════════════════════════\n\n');

% ── Per-frame rain attenuation sampler ───────────────────────────────────────
%  The paper models each frame as an independent draw from the P.618
%  complementary CDF.  We invert the P.618 power-law scaling:
%      A(p) = A001 * (p/0.01)^exponent
%  by drawing p uniformly from [p_min, p_max] and computing A(p).
%  p_min = 0.001 % → deep fade (A > A001)
%  p_max = 1.0   % → near clear sky
%  This gives physically meaningful attenuation diversity per frame.
itu.p_min = 0.001;    % % of year — deep fade end
itu.p_max = 1.0;      % % of year — near clear-sky end


%% ═══════════════════════════════════════════════════════════════════════════
%  SEC 3 — LINK BUDGET & COMPOSITE CHANNEL COEFFICIENT
%  Implements Equations (7)–(11) of the paper.
%  h_k = h̃ · b^(1/2) · sqrt(b_max)
%  where:
%    h̃      = Ap^(1/2) · exp(-j·φ)          rain fading   (Eq. 8)
%    b       = Gs · [J1(u)/2u + 36·J3(u)/u³]²  beam gain   (Eq. 10)
%    b_max   = (λ/4π)² · (1/d0²) · Gr/(kb·BW·T)  path+thermal (Eq. 9)
% ═══════════════════════════════════════════════════════════════════════════

% ── Physical constants ────────────────────────────────────────────────────────
CONST.c_light  = 3e8;                          % speed of light  [m/s]
CONST.kb       = 1.38064852e-23;               % Boltzmann const [J/K]

% ── Satellite / terminal link parameters (paper Table 4) ─────────────────────
link.freq_hz        = itu.Frequency_GHz * 1e9;
link.lambda         = CONST.c_light / link.freq_hz;
link.d0_m           = 35788e3;                 % GEO altitude    [m]
link.Gs_dBi         = 52.0;                    % satellite TX antenna gain
link.Gr_dBi         = 41.7;                    % user terminal RX antenna gain
link.NoiseBW_Hz     = 50e6;                    % noise bandwidth [Hz]
link.NoiseTemp_K    = 207;                     % clear-sky noise temperature [K]
link.Theta3dB_deg   = 0.4;                     % 3-dB beamwidth half-angle
link.BeamOffset_deg = 0.0;                     % terminal offset from beam centre

% ── Derived link quantities ───────────────────────────────────────────────────
Gr_lin   = 10^(link.Gr_dBi / 10);
Gs_lin   = 10^(link.Gs_dBi / 10);

% b_max : free-space loss + receiver noise normalisation  (Eq. 9)
b_max = (link.lambda / (4*pi))^2 ...
        * (1 / link.d0_m^2) ...
        * (Gr_lin / (CONST.kb * link.NoiseBW_Hz * link.NoiseTemp_K));

% b_gain : satellite beam pattern gain at terminal location (Eq. 10)
theta3dB_rad = deg2rad(link.Theta3dB_deg);
theta_off_rad = deg2rad(link.BeamOffset_deg);
u_val = 2.07123 * sin(theta_off_rad) / sin(theta3dB_rad);
if abs(u_val) < 1e-8
    bp = 1;    % at beam centre: pattern = 1
else
    bp = besselj(1, u_val)/(2*u_val) + 36*besselj(3, u_val)/(u_val^3);
end
b_gain = Gs_lin * abs(bp)^2;

fprintf('  Link Budget Summary\n');
fprintf('  ───────────────────────────────────────────\n');
fprintf('  Frequency     : %.0f GHz\n',   link.freq_hz/1e9);
fprintf('  GEO distance  : %.0f km\n',    link.d0_m/1e3);
fprintf('  Gs            : %.1f dBi\n',   link.Gs_dBi);
fprintf('  Gr            : %.1f dBi\n',   link.Gr_dBi);
fprintf('  Noise BW      : %.0f MHz\n',   link.NoiseBW_Hz/1e6);
fprintf('  Noise Temp    : %.0f K\n',     link.NoiseTemp_K);
fprintf('  b_max         : %.4e\n',       b_max);
fprintf('  b_gain        : %.4e\n',       b_gain);
fprintf('  ───────────────────────────────────────────\n\n');


%% ═══════════════════════════════════════════════════════════════════════════
%  SEC 4 — DVB-S2 PHY GEOMETRY  (frame structure, pilot indices)
% ═══════════════════════════════════════════════════════════════════════════

% ── Validate MODCOD selection ─────────────────────────────────────────────────
assert(all(ismember(modcod_selection, 1:13)), ...
    'Invalid MODCOD ID. Valid range is 1–13 for normal FEC frame.');
modcod_pool = modcod_selection(:)';

% ── Constants for DVB-S2 pilot structure ──────────────────────────────────────
PLHDR_LEN   = 90;      % PL header length (symbols)
SLOT_LEN    = 90;      % slot length (symbols)
PILOT_FREQ  = 16;      % pilot block every 16 slots
PILOT_SYMS  = 36;      % pilot block length (symbols)
HasPilots   = true;
syncBits    = [0;1;0;0;0;1;1;1];
pktLen      = 1496;
filterDelay = FilterSpanInSymbols / 2;

% ── Modulation order name map ─────────────────────────────────────────────────
modName = containers.Map([4 8 16 32], {'QPSK','8PSK','16APSK','32APSK'});

% ── Pre-compute geometry for each selected MODCOD ────────────────────────────
geom = containers.Map('KeyType','int32','ValueType','any');

for mc = modcod_pool
    [modOrd, cRate, cwLen] = satcom.internal.dvbs.getS2PHYParams(mc, FECframe);

    % Number of data symbols in the PLFRAME
    dLen     = cwLen / log2(modOrd);

    % Pilot block count (DVB-S2 standard: pilot every 16 slots, not at end)
    nPBlks = floor(dLen / (SLOT_LEN * PILOT_FREQ));
    if floor(dLen/(SLOT_LEN*PILOT_FREQ)) == dLen/(SLOT_LEN*PILOT_FREQ)
        nPBlks = nPBlks - 1;
    end

    pilotLen  = nPBlks * PILOT_SYMS;
    frmSz     = dLen + pilotLen + PLHDR_LEN;   % full PLFRAME length (symbols)

    % PL scrambling sequence (used to compute reference pilots)
    plScr = satcom.internal.dvbs.plScramblingIntegerSequence(0);
    cMap  = [1; 1j; -1; -1j];
    cSeq  = cMap(plScr + 1);

    % Pilot indices within PLFRAME and reference pilot values
    [~, pInd_raw] = satcom.internal.dvbs.pilotBlock(nPBlks);
    pInd_raw = pInd_raw(:)';
    pInd     = pInd_raw + PLHDR_LEN;           % offset past PL header
    refP     = (1+1j)/sqrt(2) .* cSeq(pInd_raw(:));   % known reference pilots

    g.modOrder      = modOrd;
    g.codeRate      = cRate;
    g.cwLen         = cwLen;
    g.plFrameSize   = frmSz;
    g.numPilotBlks  = nPBlks;
    g.pilotInd      = pInd;          % [Np × 1] indices into PLFRAME
    g.pilotInd_raw  = pInd_raw;
    g.refPilots     = refP;          % [Np × 1] known complex pilots
    g.Np            = length(pInd);  % total pilot symbols per frame
    g.numPkts       = [];            % filled below after wavegen init

    [cr_num, cr_den] = rat(cRate, 1e-4);
    g.label = sprintf('%s_%d_%d', modName(modOrd), cr_num, cr_den);

    geom(int32(mc)) = g;

    fprintf('  MODCOD %2d  %-10s  cwLen=%6d  frameSize=%6d  pilots=%4d\n', ...
        mc, g.label, cwLen, frmSz, g.Np);
end
fprintf('\n');

% ── Initialise wavegen once to get MinNumPackets per MODCOD ──────────────────
wavegen = dvbs2WaveformGenerator;
wavegen.StreamFormat        = 'TS';
wavegen.NumInputStreams      = 1;
wavegen.FECFrame            = FECframe;
wavegen.HasPilots            = HasPilots;
wavegen.RolloffFactor        = RolloffFactor;
wavegen.FilterSpanInSymbols  = FilterSpanInSymbols;
wavegen.SamplesPerSymbol     = SPS;

for mc = modcod_pool
    wavegen.MODCOD = mc;
    wavegen.DFL    = getDFL(mc, FECframe);
    g              = geom(int32(mc));
    g.numPkts      = wavegen.MinNumPackets(1);
    geom(int32(mc)) = g;
end

% ── Maximum frame size across all MODCODs (for zero-padding in CSV) ───────────
maxFrameSize = max(cellfun(@(mc) geom(int32(mc)).plFrameSize, num2cell(modcod_pool)));
maxNp        = max(cellfun(@(mc) geom(int32(mc)).Np,          num2cell(modcod_pool)));

% Xin has 4 × maxFrameSize columns: [Re(y), Re(p_ref), Im(y), Im(p_ref)]
Xin_cols = 4 * maxFrameSize;

fprintf('  Max PLFRAME size : %d symbols\n', maxFrameSize);
fprintf('  Max pilots (Np)  : %d\n',         maxNp);
fprintf('  Xin width        : %d values per frame\n\n', Xin_cols);


%% ═══════════════════════════════════════════════════════════════════════════
%  SEC 5 — OUTPUT FILE SETUP
% ═══════════════════════════════════════════════════════════════════════════

% ── Build output filenames ────────────────────────────────────────────────────
mc_tags = cell(1, numel(modcod_pool));
for ii = 1:numel(modcod_pool)
    mc_tags{ii} = geom(int32(modcod_pool(ii))).label;
end
mc_tag = strjoin(mc_tags, '__');

if ~exist(outDir, 'dir'), mkdir(outDir); end
matFile = fullfile(outDir, sprintf('dvbs2_%s_%s.mat', itu.channel_type, mc_tag));
csvFile = fullfile(outDir, sprintf('dvbs2_%s_%s.csv', itu.channel_type, mc_tag));

% ── Pre-allocate dataset struct ───────────────────────────────────────────────
%  Each entry stores all quantities needed to reconstruct Xin and metadata.
dataset(numFrames) = struct( ...
    'frame_idx',    [], ...  % frame counter
    'modcod',       [], ...  % MODCOD ID (integer)
    'modcod_label', [], ...  % e.g. 'QPSK_1_2'
    'channel_type', [], ...  % 'tropical'
    'esno_nom_dB',  [], ...  % nominal Es/No drawn from esno_pool (before rain)
    'rainAtt_dB',   [], ...  % ITU-R P.618 rain attenuation this frame  [dB]
    'p_exceedance', [], ...  % exceedance probability p used for this frame [%]
    'snr_eff_dB',   [], ...  % effective Es/No = esno_nom - rainAtt        [dB]
    'nvar',         [], ...  % noise variance estimated from received pilots
    'Ap_linear',    [], ...  % rain amplitude attenuation factor (linear, >1 = loss)
    'phi_rad',      [], ...  % random phase rotation applied [rad]
    'b_gain',       [], ...  % beam pattern gain (scalar, same every frame)
    'b_max',        [], ...  % free-space + thermal normalisation (scalar)
    'htrue_real',   [], ...  % real(H_true) — composite channel coefficient
    'htrue_imag',   [], ...  % imag(H_true)
    'Xin_real_y',   [], ...  % real(y_full) — [plFrameSize × 1]
    'Xin_imag_y',   [], ...  % imag(y_full) — [plFrameSize × 1]
    'Xin_real_p',   [], ...  % real(p_ref)  — [plFrameSize × 1] zero-padded
    'Xin_imag_p',   []);     % imag(p_ref)  — [plFrameSize × 1] zero-padded

% ── Open CSV and write header ─────────────────────────────────────────────────
if SAVE_CSV
    fid = fopen(csvFile, 'w');
    if fid == -1, error('Cannot open CSV file: %s', csvFile); end

    % Metadata columns
    hdr = ['frame_idx,modcod,modcod_label,channel_type,' ...
           'esno_nom_dB,rainAtt_dB,p_exceedance,snr_eff_dB,' ...
           'nvar,Ap_linear,phi_rad,b_gain,b_max,' ...
           'htrue_real,htrue_imag,'];

    % Xin columns: Re(y_1)…Re(y_N), Re(p_1)…Re(p_N), Im(y), Im(p)
    for i = 1:maxFrameSize, hdr = [hdr sprintf('Xin_rey_%d,',  i)]; end %#ok
    for i = 1:maxFrameSize, hdr = [hdr sprintf('Xin_rep_%d,',  i)]; end %#ok
    for i = 1:maxFrameSize, hdr = [hdr sprintf('Xin_imy_%d,',  i)]; end %#ok
    for i = 1:maxFrameSize, hdr = [hdr sprintf('Xin_imp_%d,',  i)]; end %#ok
    hdr(end) = [];   % remove trailing comma

    fprintf(fid, '%s\n', hdr);
    fprintf('  CSV header written  (%d metadata + %d Xin columns)\n\n', ...
        15, Xin_cols);
end


%% ═══════════════════════════════════════════════════════════════════════════
%  MAIN GENERATION LOOP
% ═══════════════════════════════════════════════════════════════════════════

fprintf('Generating %d frames …\n\n', numFrames);
fprintf('%-7s %-12s %-14s %-13s %-13s %-12s %-10s\n', ...
    'Frame','MODCOD','EsNo_nom(dB)','RainAtt(dB)','EsNo_eff(dB)','p_exceed(%)','nvar');
fprintf('%s\n', repmat('─', 1, 83));

for frameIdx = 1 : numFrames

    %% ─────────────────────────────────────────────────────────────────────
    %  SEC 6 — TRANSMITTER
    %  Randomly pick MODCOD and Es/No for this frame, then build PLFRAME.
    %% ─────────────────────────────────────────────────────────────────────

    mc      = modcod_pool(randi(length(modcod_pool)));
    esno_dB = esno_pool(randi(length(esno_pool)));
    g       = geom(int32(mc));

    % Generate transport stream bits (PN sequence, same as original)
    pn = comm.PNSequence('Polynomial',        'x9+x5+1', ...
                         'InitialConditions',  [zeros(1,8) 1], ...
                         'VariableSizeOutput', true, ...
                         'MaximumOutputSize',  [pktLen * g.numPkts, 1]);
    reset(pn);
    rawBits   = pn(pktLen * g.numPkts);
    txRawPkts = reshape(rawBits, pktLen, g.numPkts);
    txPkts    = [repmat(syncBits, 1, g.numPkts); txRawPkts];
    inputBits = txPkts(:);

    release(wavegen);
    wavegen.MODCOD = mc;
    wavegen.DFL    = getDFL(mc, FECframe);
    txWaveform     = [wavegen({inputBits}); flushFilter(wavegen)];

    clearvars rawBits txRawPkts txPkts inputBits pn;


    %% ─────────────────────────────────────────────────────────────────────
    %  SEC 7 — CHANNEL
    %
    %  Step 1 : Draw rain attenuation from ITU-R P.618 distribution
    %           by sampling exceedance probability p uniformly in [p_min, p_max]
    %           and computing A(p) via the P.618 power-law scaling.
    %
    %  Step 2 : Build composite channel coefficient h_k (Eq. 11):
    %           h_k = h̃ · sqrt(b_gain) · sqrt(b_max)
    %           h̃   = (1/sqrt(Ap)) · exp(-j·φ)   (Eq. 8)
    %           Ap  = 10^(rainAtt_dB / 20)        (linear amplitude factor)
    %
    %  Step 3 : Effective Es/No = esno_nom - rainAtt_dB
    %           Noise variance σ²_n set from effective Es/No.
    %           AWGN added at oversampled (SPS=2) waveform level.
    %% ─────────────────────────────────────────────────────────────────────

    % -- 7a. ITU-R P.618 rain attenuation this frame -------------------------
    %  Draw p uniformly in [p_min, p_max], then compute A(p).
    %  p is in percent-of-year.  Lower p → rarer event → heavier rain.
    p_draw      = itu.p_min + (itu.p_max - itu.p_min) * rand();
    rainAtt_dB  = calcRainAttenuationP618(itu, p_draw);
    rainAtt_dB  = max(0, rainAtt_dB);   % guard against numerical negatives

    % -- 7b. Effective received Es/No ----------------------------------------
    snr_eff_dB = esno_dB - rainAtt_dB;

    % SNR seen by AWGN block must account for oversampling (SPS):
    %   Es/No_awgn = Es/No_eff - 10·log10(SPS)
    snr_awgn = snr_eff_dB - 10*log10(SPS);

    % -- 7c. Composite channel coefficient h_k (Eq. 7–11) -------------------
    Ap            = 10^(rainAtt_dB / 20);         % amplitude attenuation (>1)
    phi           = 2 * pi * rand();               % uniform phase rotation [rad]
    h_tilde       = (1 / sqrt(Ap)) * exp(-1j * phi);  % rain fading (Eq. 8)
    H_true        = h_tilde * sqrt(b_gain) * sqrt(b_max);  % full h_k (Eq. 11)

    % -- 7d. Apply channel to waveform ----------------------------------------
    attWaveform = H_true * txWaveform;
    rxNoisy     = awgn(attWaveform, snr_awgn, 'measured');

    clearvars txWaveform attWaveform;


    %% ─────────────────────────────────────────────────────────────────────
    %  SEC 8 — RECEIVER
    %  Matched filter → symbol timing → extract PLFRAME symbols.
    %% ─────────────────────────────────────────────────────────────────────

    % -- 8a. Matched filter (root-raised cosine receive filter) ---------------
    rxFilter = comm.RaisedCosineReceiveFilter( ...
        'RolloffFactor',         RolloffFactor, ...
        'FilterSpanInSymbols',   FilterSpanInSymbols, ...
        'InputSamplesPerSymbol', SPS, ...
        'DecimationFactor',      SPS);
    rxSym_raw = rxFilter(rxNoisy);
    release(rxFilter);
    clearvars rxNoisy;

    % -- 8b. Remove filter transient (symmetric truncation) -------------------
    rxSymbols = rxSym_raw(filterDelay + 1 : end - filterDelay);
    clearvars rxSym_raw;

    % -- 8c. Guard: enough symbols? -------------------------------------------
    if length(rxSymbols) < g.plFrameSize
        warning('Frame %d: too few symbols after filtering — skipping.', frameIdx);
        clearvars rxSymbols;
        continue;
    end

    % -- 8d. Extract first complete PLFRAME -----------------------------------
    rxFrame = rxSymbols(1 : g.plFrameSize);   % [plFrameSize × 1] complex
    clearvars rxSymbols;

    % -- 8e. Pilot-based noise variance estimate ------------------------------
    %  Uses received pilots vs known reference pilots.
    %  nvar = E[|rx_pilot|²] - |E[rx_pilot · conj(ref_pilot)]|²
    nvar = pilotNoiseVarEstimate(rxFrame, g.pilotInd, g.refPilots);


    %% ─────────────────────────────────────────────────────────────────────
    %  SEC 9 — Xin ASSEMBLY & SAVING
    %
    %  Paper Eq. 33 (LS part omitted as requested):
    %      Xin = [ real(y)^T,  real(p_ref_padded)^T,
    %              imag(y)^T,  imag(p_ref_padded)^T ]^T
    %
    %  y         = full received PLFRAME  [plFrameSize × 1]
    %  p_ref_pad = reference pilots placed at pilot indices,
    %              zeros at data positions [plFrameSize × 1]
    %% ─────────────────────────────────────────────────────────────────────

    % -- 9a. Build p_ref_padded (frame-length vector, pilots at their positions)
    p_ref_padded = zeros(g.plFrameSize, 1);
    p_ref_padded(g.pilotInd(:)) = g.refPilots(:);

    % -- 9b. Xin components (as paper defines, zero-padded to maxFrameSize) ---
    Xin_rey = zeros(maxFrameSize, 1);
    Xin_imy = zeros(maxFrameSize, 1);
    Xin_rep = zeros(maxFrameSize, 1);
    Xin_imp = zeros(maxFrameSize, 1);

    Xin_rey(1:g.plFrameSize) = real(rxFrame);
    Xin_imy(1:g.plFrameSize) = imag(rxFrame);
    Xin_rep(1:g.plFrameSize) = real(p_ref_padded);
    Xin_imp(1:g.plFrameSize) = imag(p_ref_padded);

    % -- 9c. Store in dataset struct ------------------------------------------
    dataset(frameIdx).frame_idx    = frameIdx;
    dataset(frameIdx).modcod       = mc;
    dataset(frameIdx).modcod_label = g.label;
    dataset(frameIdx).channel_type = itu.channel_type;
    dataset(frameIdx).esno_nom_dB  = esno_dB;
    dataset(frameIdx).rainAtt_dB   = rainAtt_dB;
    dataset(frameIdx).p_exceedance = p_draw;
    dataset(frameIdx).snr_eff_dB   = snr_eff_dB;
    dataset(frameIdx).nvar         = nvar;
    dataset(frameIdx).Ap_linear    = Ap;
    dataset(frameIdx).phi_rad      = phi;
    dataset(frameIdx).b_gain       = b_gain;
    dataset(frameIdx).b_max        = b_max;
    dataset(frameIdx).htrue_real   = real(H_true);
    dataset(frameIdx).htrue_imag   = imag(H_true);
    dataset(frameIdx).Xin_real_y   = Xin_rey;
    dataset(frameIdx).Xin_imag_y   = Xin_imy;
    dataset(frameIdx).Xin_real_p   = Xin_rep;
    dataset(frameIdx).Xin_imag_p   = Xin_imp;

    % -- 9d. Write CSV row ----------------------------------------------------
    if SAVE_CSV
        metaVec = [frameIdx, mc, NaN, NaN, ...   % label fields handled as strings
                   esno_dB, rainAtt_dB, p_draw, snr_eff_dB, ...
                   nvar, Ap, phi, b_gain, b_max, ...
                   real(H_true), imag(H_true)];

        % Build numeric part of row (metadata + Xin flattened)
        XinFlat = [Xin_rey; Xin_rep; Xin_imy; Xin_imp]';   % 1 × (4×maxFrameSize)

        % Write metadata string columns first, then numeric
        fprintf(fid, '%d,%d,%s,%s,', frameIdx, mc, g.label, itu.channel_type);
        fprintf(fid, '%.6g,%.6g,%.8g,%.6g,', esno_dB, rainAtt_dB, p_draw, snr_eff_dB);
        fprintf(fid, '%.8g,%.8g,%.8g,%.8g,%.8g,', nvar, Ap, phi, b_gain, b_max);
        fprintf(fid, '%.8g,%.8g,', real(H_true), imag(H_true));

        % Write Xin values — use compact format to keep file size manageable
        fprintf(fid, '%.6g,', XinFlat(1:end-1));
        fprintf(fid, '%.6g\n', XinFlat(end));
    end

    % -- 9e. Console progress -------------------------------------------------
    fprintf('%-7d %-12s %-14.2f %-13.4f %-13.4f %-12.6f %-10.4e\n', ...
        frameIdx, g.label, esno_dB, rainAtt_dB, snr_eff_dB, p_draw, nvar);

    % -- 9f. Cleanup ----------------------------------------------------------
    clearvars rxFrame p_ref_padded Xin_rey Xin_imy Xin_rep Xin_imp XinFlat ...
              nvar H_true h_tilde Ap phi rainAtt_dB p_draw snr_eff_dB ...
              snr_awgn g mc esno_dB metaVec;

    % -- 9g. MAT checkpoint ---------------------------------------------------
    if mod(frameIdx, saveEvery) == 0 || frameIdx == numFrames
        save(matFile, 'dataset', 'itu', 'link', 'CONST', '-v7.3');
        fprintf('  → Checkpoint %d/%d saved → %s\n', frameIdx, numFrames, matFile);
    end

end % ── end main loop ────────────────────────────────────────────────────────

if SAVE_CSV, fclose(fid); end

fprintf('\n══════════════════════════════════════════════════════\n');
fprintf('  Dataset generation complete\n');
fprintf('  Frames      : %d\n',  numFrames);
fprintf('  MAT file    → %s\n',  matFile);
if SAVE_CSV
    fprintf('  CSV file    → %s\n',  csvFile);
end
fprintf('══════════════════════════════════════════════════════\n');


%% ══════════════════════════════════════════════════════════════════════════
%  HELPER FUNCTIONS
%  Each function is self-contained with a clear description.
%% ══════════════════════════════════════════════════════════════════════════


function dfl = getDFL(modCod, fecFrame)
% GETDFL  Data field length (DFL) for a given MODCOD + FEC frame type.
%   Used to set wavegen.DFL before generating the waveform.
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


function nVarEst = pilotNoiseVarEstimate(rxFrame, pilotInd, refPilots)
% PILOTNOISEVARESTIMATE  Estimate noise variance from received pilot symbols.
%
%   Method: for a received pilot r = h·p + n where h≈constant and p is known,
%     signal power  = |E[r · conj(p)]|²  =  |h|²·|p|²
%     total power   = E[|r|²]            =  |h|²·|p|² + σ²_n
%     noise var est = total power − signal power
%
%   Inputs:
%     rxFrame   — full received PLFRAME  [N × 1] complex
%     pilotInd  — indices of pilot symbols in rxFrame
%     refPilots — known reference pilot symbols [Np × 1] complex
%
%   Output:
%     nVarEst  — estimated noise variance (scalar, ≥ 0)

rxPilots  = rxFrame(pilotInd(:));
totalPow  = mean(abs(rxPilots).^2);
sigPow    = abs(mean(rxPilots .* conj(refPilots(:)))).^2;
nVarEst   = max(0, totalPow - sigPow);
end


function Ap_dB = calcRainAttenuationP618(chanParams, p)
% CALCRAINSATTENUATIONP618  ITU-R P.618 rain attenuation at exceedance p (%).
%
%   Implements the full P.618-13 prediction method:
%     Step 1 : Compute slant path length Ls  (ITU-R P.839 rain height)
%     Step 2 : Specific attenuation γ_R      (ITU-R P.838 k, α coefficients)
%     Step 3 : Horizontal / vertical reduction factors r001, v001
%     Step 4 : Effective path length LE → A_001
%     Step 5 : Scale A_001 to arbitrary p via power-law
%
%   Inputs:
%     chanParams — struct with fields:
%       .Frequency_GHz, .Latitude_deg, .Altitude_m, .ElevationAngle_deg,
%       .Polarization ('H'|'V'), .R001_mmph, .IsothermHeight_km
%     p          — exceedance probability [% of year], e.g. 0.01
%
%   Output:
%     Ap_dB      — rain attenuation exceeded p% of the year [dB]

f_GHz     = chanParams.Frequency_GHz;
theta_deg = chanParams.ElevationAngle_deg;
lat_deg   = chanParams.Latitude_deg;
hs_km     = chanParams.Altitude_m / 1000;
R001      = chanParams.R001_mmph;
h0_km     = chanParams.IsothermHeight_km;

% Rain height (ITU-R P.839)
hR_km = h0_km + 0.36;
Re_km = 8500;   % effective Earth radius [km]

if R001 <= 0 || hR_km <= hs_km
    Ap_dB = 0; return;
end

% Slant path length Ls [km]
if theta_deg >= 5
    Ls = (hR_km - hs_km) / sind(theta_deg);
else
    Ls = 2*(hR_km - hs_km) / ...
         (sqrt(sind(theta_deg)^2 + 2*(hR_km-hs_km)/Re_km) + sind(theta_deg));
end
LG = Ls * cosd(theta_deg);   % horizontal projection [km]

% Specific attenuation coefficients k, α  (ITU-R P.838)
[k, alpha] = ituP838_k_alpha(f_GHz, theta_deg, chanParams.Polarization);
gamma_R    = k * (R001^alpha);   % specific attenuation [dB/km] at R001

% Horizontal reduction factor r001
r001 = 1 / (1 + 0.78*sqrt((LG*gamma_R)/f_GHz) - 0.38*(1 - exp(-2*LG)));

% Effective slant path length
zeta = atand((hR_km - hs_km) / (LG*r001));
if zeta > theta_deg
    LR = (LG*r001) / cosd(theta_deg);
else
    LR = (hR_km - hs_km) / sind(theta_deg);
end

% Vertical adjustment factor v001
chi  = max(0, 36 - abs(lat_deg));
v001 = 1 / (1 + sqrt(sind(theta_deg)) * ...
       (31*(1 - exp(-theta_deg/(1+chi))) * sqrt(LR*gamma_R)/(f_GHz^2) - 0.45));

LE   = LR * v001;       % effective path length [km]
A001 = gamma_R * LE;    % attenuation at 0.01% exceedance [dB]

% Return A001 directly if p == 0.01%
if abs(p - 0.01) < 1e-12
    Ap_dB = A001; return;
end

% Scale to arbitrary p using P.618 power-law
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

exponent = -(0.655 + 0.033*log(p) - 0.045*log(A001) ...
             - beta*(1-p)*sind(theta_deg));
Ap_dB    = A001 * (p/0.01)^exponent;
end


function [k, alpha] = ituP838_k_alpha(f_GHz, theta_deg, polarization)
% ITUP838_K_ALPHA  Specific attenuation model coefficients (ITU-R P.838-3).
%
%   Returns the frequency- and polarisation-dependent coefficients k and α
%   used in:  γ_R = k · R^α  [dB/km]
%
%   Inputs:
%     f_GHz       — frequency [GHz]
%     theta_deg   — path elevation angle [degrees]
%     polarization — 'H' (horizontal) | 'V' (vertical)
%
%   Output:
%     k, alpha    — specific attenuation coefficients (scalars)

x = log10(f_GHz);

% Regression coefficients for k_H and k_V
akH = [-5.33980, -0.35351, -0.23789, -0.94158];
bkH = [-0.10008,  1.26970,  0.86036,  0.64552];
ckH = [ 1.13098,  0.45400,  0.15354,  0.16817];

akV = [-3.80595, -3.44965, -0.39902,  0.50167];
bkV = [ 0.56934, -0.22911,  0.73042,  1.07319];
ckV = [ 0.81061,  0.51059,  0.11899,  0.27195];

logkH = sum(akH .* exp(-((x-bkH)./ckH).^2)) - 0.18961*x + 0.71147;
logkV = sum(akV .* exp(-((x-bkV)./ckV).^2)) - 0.16398*x + 0.63297;
kH = 10^logkH;
kV = 10^logkV;

% Regression coefficients for α_H and α_V
aaH = [-0.14318,  0.29591,  0.32177, -5.37610, 16.1721];
baH = [ 1.82442,  0.77564,  0.63773, -0.96230, -3.29980];
caH = [-0.55187,  0.19822,  0.13164,  1.47828,  3.43990];

aaV = [-0.07771,  0.56727, -0.20238, -48.2991, 48.5833];
baV = [ 2.33840,  0.95545,  1.14520,  0.791669, 0.791459];
caV = [-0.76284,  0.54039,  0.26809,  0.116226, 0.116479];

alphaH = sum(aaH .* exp(-((x-baH)./caH).^2)) + 0.67849*x - 1.95537;
alphaV = sum(aaV .* exp(-((x-baV)./caV).^2)) - 0.053739*x + 0.83433;

% Combined k and α for arbitrary polarisation and elevation (P.838 Eq. 4-5)
pol = upper(string(polarization));
if     pol == "H", tau_deg = 0;
elseif pol == "V", tau_deg = 90;
else,              tau_deg = 45;
end

k     = (kH + kV + (kH-kV)*cosd(theta_deg)^2*cosd(2*tau_deg)) / 2;
alpha = (kH*alphaH + kV*alphaV + ...
        (kH*alphaH - kV*alphaV)*cosd(theta_deg)^2*cosd(2*tau_deg)) / (2*k);
end