% ========================= DVB-S2 Per-MODCOD Dataset Generator =========================
%
% OUTPUT per row:
%   [ metadata columns ] + [ pilot_re_1…MAX_NP ] + [ pilot_im_1…MAX_NP ] + [ pilot_mask_1…MAX_NP ]
%
% HOW PILOT PADDING WORKS:
%   Each MODCOD has a fixed Np (pilot symbol count).  Because every file
%   contains only ONE MODCOD, MAX_NP equals that MODCOD's own Np — no
%   cross-MODCOD padding is needed within a single file.  The mask column
%   is kept for consistency with the multi-MODCOD generator (all ones).
%
% CHANNEL:  ITU-R P.618 rain-fading + AWGN (same model as dataset_gen.m,
%           Eq. 7-11).  Es/No is swept deterministically (samples_per_esno
%           frames per level), in the same spirit as ASK2_data_generation.m.


clc; clearvars; close all;

%% 
% USER KNOBS  ◄ EDIT ONLY THIS SECTION FOR NORMAL USE ►

modcod_selection = [2, 9, 15, 16, 19, 22, 24, 26];    % ← replace with your 8 IDs

% ── Es/No sweep ──────────────────────────────────────────────────────────────
%   Deterministic: samples_per_esno frames are generated for EACH Es/No level.
esno_values     = -2:2:18;          % dB  (same style as ASK2 snr_values)
samples_per_esno = 1;            % frames per Es/No level per MODCOD

% ── Output ───────────────────────────────────────────────────────────────────
outDir = '../dataset_output';       % folder (created if absent)
%   Each MODCOD writes its own file:  <outDir>/<modcod_label>.csv
%   e.g.  QPSK_1_2.csv,  8PSK_3_4.csv,  …

% ── FEC frame type ───────────────────────────────────────────────────────────
FECframe = 'normal';                % 'normal' (64 800 b) | 'short' (16 200 b)

% ── Waveform / filter ────────────────────────────────────────────────────────
RolloffFactor       = 0.20;
FilterSpanInSymbols = 10;
SPS                 = 2;            % samples per symbol

% ── Random seed  (0 = random, else fixed for reproducibility) ───────────────
RNG_SEED = 0;

% ── Debug print  (one line per first sample of each MODCOD / Es/No pair) ────
enableDebugPrint = true;


%% ═══════════════════════════════════════════════════════════════════════════
%  SEC 2 — CHANNEL PLAYGROUND  (ITU-R P.618 rain-fading, same as dataset_gen)
%% ═══════════════════════════════════════════════════════════════════════════

itu.channel_type        = 'tropical';
itu.site                = 'Penang, Malaysia';
itu.Frequency_GHz       = 18;
itu.Latitude_deg        = 23.999324;
itu.Longitude_deg       = 90.389147;
itu.Altitude_m          = 57;
itu.ElevationAngle_deg  = 45;
itu.Polarization        = 'V';
itu.R001_mmph           = 130;
itu.IsothermHeight_km   = 4.5;
itu.p_min               = 0.01;    % deep-fade end
itu.p_max               = 1.0;     % near clear-sky end

link.Gs_dBi         = 52.0;
link.Gr_dBi         = 41.7;
link.NoiseBW_Hz     = 50e6;
link.NoiseTemp_K    = 207;
link.Theta3dB_deg   = 0.4;
link.BeamOffset_deg = 0.0;
link.d0_m           = 35788e3;

CONST.c_light = 3e8;
CONST.kb      = 1.38064852e-23;





if RNG_SEED ~= 0, rng(RNG_SEED); end

link.freq_hz = itu.Frequency_GHz * 1e9;
link.lambda  = CONST.c_light / link.freq_hz;

%% ── Link-budget derived quantities (Eq. 9-10) ───────────────────────────────
Gr_lin = 10^(link.Gr_dBi / 10);
Gs_lin = 10^(link.Gs_dBi / 10);

b_max = (link.lambda / (4*pi))^2 * (1 / link.d0_m^2) ...
        * (Gr_lin / (CONST.kb * link.NoiseBW_Hz * link.NoiseTemp_K));

theta3dB_rad  = deg2rad(link.Theta3dB_deg);
theta_off_rad = deg2rad(link.BeamOffset_deg);
u_val = 2.07123 * sin(theta_off_rad) / sin(theta3dB_rad);
if abs(u_val) < 1e-8
    bp = 1;
else
    bp = besselj(1, u_val)/(2*u_val) + 36*besselj(3, u_val)/(u_val^3);
end
b_gain = Gs_lin * abs(bp)^2;

A001_dB = calcRainAttenuationP618(itu, 0.01);

fprintf('═══════════════════════════════════════════════════════\n');
fprintf('  ITU-R P.618  |  %s  |  A_001 = %.2f dB\n', itu.site, A001_dB);
fprintf('  b_max = %.4e   b_gain = %.4e\n', b_max, b_gain);
fprintf('═══════════════════════════════════════════════════════\n\n');


%% ═══════════════════════════════════════════════════════════════════════════
%  SEC 3 — DVB-S2 PHY GEOMETRY (one entry per selected MODCOD)
%% ═══════════════════════════════════════════════════════════════════════════

assert(all(ismember(modcod_selection, 1:28)), ...
    'Invalid MODCOD ID(s). Valid range is 1–28.');

PLHDR_LEN  = 90;
SLOT_LEN   = 90;
PILOT_FREQ = 16;
PILOT_SYMS = 36;
HasPilots  = true;
syncBits   = [0;1;0;0;0;1;1;1];
pktLen     = 1496;
filterDelay = FilterSpanInSymbols / 2;

modName = containers.Map([4 8 16 32], {'QPSK','8PSK','16APSK','32APSK'});
geom    = containers.Map('KeyType','int32','ValueType','any');

fprintf('  MODCOD geometry:\n');
for mc = modcod_selection(:)'
    [modOrd, cRate, cwLen] = satcom.internal.dvbs.getS2PHYParams(mc, FECframe);

    dLen   = cwLen / log2(modOrd);
    nPBlks = floor(dLen / (SLOT_LEN * PILOT_FREQ));
    if floor(dLen/(SLOT_LEN*PILOT_FREQ)) == dLen/(SLOT_LEN*PILOT_FREQ)
        nPBlks = nPBlks - 1;
    end

    pilotLen = nPBlks * PILOT_SYMS;
    frmSz    = dLen + pilotLen + PLHDR_LEN;

    plScr = satcom.internal.dvbs.plScramblingIntegerSequence(0);
    cMap  = [1; 1j; -1; -1j];
    cSeq  = cMap(plScr + 1);

    [~, pInd_raw] = satcom.internal.dvbs.pilotBlock(nPBlks);
    pInd_raw = pInd_raw(:)';
    pInd     = pInd_raw + PLHDR_LEN;
    refP     = (1+1j)/sqrt(2) .* cSeq(pInd_raw(:));

    [cr_num, cr_den] = rat(cRate, 1e-4);

    g.modOrder     = modOrd;
    g.codeRate     = cRate;
    g.cwLen        = cwLen;
    g.plFrameSize  = frmSz;
    g.numPilotBlks = nPBlks;
    g.pilotInd     = pInd(:);
    g.pilotInd_raw = pInd_raw(:);
    g.refPilots    = refP(:);
    g.Np           = length(pInd);
    g.numPkts      = [];
    g.cr_num       = cr_num;
    g.cr_den       = cr_den;
    g.label        = sprintf('%s_%d_%d', modName(modOrd), cr_num, cr_den);

    geom(int32(mc)) = g;

    fprintf('    MODCOD %2d  %-12s  Np=%4d  frameSize=%6d\n', ...
        mc, g.label, g.Np, frmSz);
end
fprintf('\n');

%% ── Initialise wavegen once; query MinNumPackets per MODCOD ─────────────────
wavegen = dvbs2WaveformGenerator;
wavegen.StreamFormat        = 'TS';
wavegen.NumInputStreams      = 1;
wavegen.FECFrame             = FECframe;
wavegen.HasPilots            = HasPilots;
wavegen.RolloffFactor        = RolloffFactor;
wavegen.FilterSpanInSymbols  = FilterSpanInSymbols;
wavegen.SamplesPerSymbol     = SPS;

for mc = modcod_selection(:)'
    wavegen.MODCOD = mc;
    wavegen.DFL    = getDFL(mc, FECframe);
    g              = geom(int32(mc));
    g.numPkts      = wavegen.MinNumPackets(1);
    geom(int32(mc)) = g;
end

%% ── Save reference pilot lookup (once, shared across all files) ─────────────
if ~exist(outDir, 'dir'), mkdir(outDir); end

refPilotLookup = struct();
for mc = modcod_selection(:)'
    g  = geom(int32(mc));
    fn = sprintf('modcod_%d', mc);
    refPilotLookup.(fn).label     = g.label;
    refPilotLookup.(fn).Np        = g.Np;
    refPilotLookup.(fn).refPilots = g.refPilots;
    refPilotLookup.(fn).pilotInd  = g.pilotInd;
end
refPilotFile = fullfile(outDir, 'ref_pilots_lookup.mat');
save(refPilotFile, 'refPilotLookup', 'modcod_selection', 'FECframe');
fprintf('  Reference pilot lookup → %s\n\n', refPilotFile);


%% ═══════════════════════════════════════════════════════════════════════════
%  SEC 4 — PER-MODCOD GENERATION LOOP
%  Outer loop : MODCOD  → one CSV file per MODCOD
%  Middle loop: Es/No level (deterministic sweep, like ASK2 snr_values loop)
%  Inner loop : sample index 1..samples_per_esno
%% ═══════════════════════════════════════════════════════════════════════════

totalFrames = numel(modcod_selection) * numel(esno_values) * samples_per_esno;
fprintf('  Total frames to generate: %d\n', totalFrames);
fprintf('  (= %d MODCODs × %d Es/No levels × %d samples/level)\n\n', ...
    numel(modcod_selection), numel(esno_values), samples_per_esno);

for mc = modcod_selection(:)'

    g      = geom(int32(mc));
    MAX_NP = g.Np;          % single-MODCOD file → no cross-MODCOD padding needed

    csvFile = fullfile(outDir, sprintf('%s.csv', g.label));
    fid     = fopen(csvFile, 'w');
    if fid == -1
        error('Cannot open CSV: %s', csvFile);
    end

    % ── CSV header ────────────────────────────────────────────────────────────
    metaFields = { ...
        'sample_id', 'modcod', 'modcod_label', 'channel_type', ...
        'modOrder', 'codeRate_num', 'codeRate_den', 'plFrameSize', ...
        'Np_actual', 'numPilotBlks', ...
        'esno_nom_dB', 'rainAtt_dB', 'p_exceedance', 'snr_eff_dB', 'snr_awgn_dB', ...
        'nvar_pilot', 'pilotScale', ...
        'Ap_linear', 'phi_rad', 'b_gain', 'b_max', ...
        'htrue_real', 'htrue_imag', 'htrue_mag', 'htrue_phase_rad'};

    hdr = strjoin(metaFields, ',');
    for i = 1:MAX_NP, hdr = [hdr sprintf(',pilot_re_%d',   i)]; end %#ok<AGROW>
    for i = 1:MAX_NP, hdr = [hdr sprintf(',pilot_im_%d',   i)]; end %#ok<AGROW>
    for i = 1:MAX_NP, hdr = [hdr sprintf(',pilot_mask_%d', i)]; end %#ok<AGROW>
    fprintf(fid, '%s\n', hdr);

    % ── Phase random-walk initialisation (slow GEO phase drift) ──────────────
    phi       = 2 * pi * rand();
    sigma_phi = 0.02;               % rad / frame

    sample_id = 1;
    numMetaCols = numel(metaFields);
    PILOT_COLS  = 3 * MAX_NP;

    fprintf('══════════════════════════════════════════════════════\n');
    fprintf('  MODCOD %2d  %-12s  Np=%d  → %s\n', mc, g.label, MAX_NP, csvFile);
    fprintf('  %-8s %-7s %-13s %-12s %-12s %-12s\n', ...
        'SampleID','EsNo','RainAtt(dB)','EsNo_eff','p_exceed','nvar');
    fprintf('%s\n', repmat('─', 1, 70));

    for esno_dB = esno_values

        for s = 1 : samples_per_esno

            %% ── SEC 5: TRANSMITTER ──────────────────────────────────────────

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


            %% ── SEC 6: CHANNEL  (ITU-R P.618 rain fading + AWGN) ────────────

            % 6a. Log-uniform rain exceedance sample (better fade coverage)
            logp       = log10(itu.p_min) + rand()*(log10(itu.p_max)-log10(itu.p_min));
            p_draw     = 10^logp;
            rainAtt_dB = calcRainAttenuationP618(itu, p_draw);
            rainAtt_dB = max(0, min(rainAtt_dB, 40));   % clip 0..40 dB

            % 6b. Effective Es/No: rain decoupled from noise (see dataset_gen)
            snr_eff_dB = esno_dB;
            snr_awgn   = snr_eff_dB - 10*log10(SPS);

            % 6c. Channel coefficient H_true (Eq. 7-11)
            Ap        = 10^(rainAtt_dB / 20);
            delta_phi = sigma_phi * randn();
            phi       = phi + delta_phi;
            h_tilde   = (1 / sqrt(Ap)) * exp(-1j * phi);
            H_true    = h_tilde * sqrt(b_gain) * sqrt(b_max);

            % 6d. Apply channel + AWGN
            attWaveform = H_true * txWaveform;
            rxPower     = mean(abs(attWaveform).^2);
            noiseVar    = rxPower / (10^(snr_awgn/10));
            noise       = sqrt(noiseVar/2) * ...
                          (randn(size(attWaveform)) + 1j*randn(size(attWaveform)));
            rxNoisy     = attWaveform + noise;

            clearvars txWaveform attWaveform noise;


            %% ── SEC 7: RECEIVER (matched filter + frame extraction) ──────────

            rxFilter = comm.RaisedCosineReceiveFilter( ...
                'RolloffFactor',         RolloffFactor, ...
                'FilterSpanInSymbols',   FilterSpanInSymbols, ...
                'InputSamplesPerSymbol', SPS, ...
                'DecimationFactor',      SPS);
            rxSym_raw = rxFilter(rxNoisy);
            release(rxFilter);
            clearvars rxNoisy;

            rxSymbols = rxSym_raw(filterDelay + 1 : end - filterDelay);
            clearvars rxSym_raw;

            if length(rxSymbols) < g.plFrameSize
                warning('MODCOD %d | EsNo %.1f dB | sample %d: too few symbols — skipping.', ...
                    mc, esno_dB, s);
                clearvars rxSymbols;
                continue;
            end

            rxFrame = rxSymbols(1 : g.plFrameSize);
            clearvars rxSymbols;


            %% ── SEC 8: PILOT EXTRACTION & PADDING ───────────────────────────

            rxPilots = rxFrame(g.pilotInd);     % [Np × 1] complex

            % Normalize pilots for NN convergence (recommendation #6)
            pilotScale = max(abs(rxPilots));
            if pilotScale > 0
                rxPilots = rxPilots / pilotScale;
            end

            % Zero-pad to MAX_NP (= Np for single-MODCOD file; mask all ones)
            pilot_real = zeros(MAX_NP, 1);
            pilot_imag = zeros(MAX_NP, 1);
            pilot_mask = zeros(MAX_NP, 1);
            pilot_real(1:g.Np) = real(rxPilots);
            pilot_imag(1:g.Np) = imag(rxPilots);
            pilot_mask(1:g.Np) = 1;

            % Pilot-based noise variance estimate (metadata / diagnostic)
            nvar_pilot = pilotNoiseVarEstimate(rxFrame, g.pilotInd, g.refPilots);

            clearvars rxFrame rxPilots;


            %% ── SEC 9: ROW ASSEMBLY & CSV WRITE ─────────────────────────────

            % Metadata
            fprintf(fid, '%d,%d,%s,%s,%d,%d,%d,%d,%d,%d,', ...
                sample_id, mc, g.label, itu.channel_type, ...
                g.modOrder, g.cr_num, g.cr_den, g.plFrameSize, g.Np, g.numPilotBlks);
            fprintf(fid, '%.6g,%.6g,%.8g,%.6g,%.6g,', ...
                esno_dB, rainAtt_dB, p_draw, snr_eff_dB, snr_awgn);
            fprintf(fid, '%.8g,%.8g,', nvar_pilot, pilotScale);
            fprintf(fid, '%.8g,%.8g,%.8g,%.8g,', Ap, phi, b_gain, b_max);
            fprintf(fid, '%.8g,%.8g,%.8g,%.8g', ...
                real(H_true), imag(H_true), abs(H_true), angle(H_true));

            % Pilot feature columns
            fprintf(fid, ',%.6g', pilot_real);
            fprintf(fid, ',%.6g', pilot_imag);
            fprintf(fid, ',%d',   pilot_mask);
            fprintf(fid, '\n');

            % Debug print — first sample of each (MODCOD, Es/No) pair
            if enableDebugPrint && s == 1
                fprintf('  %-8d %-7.1f %-13.4f %-12.4f %-12.6f %-12.4e\n', ...
                    sample_id, esno_dB, rainAtt_dB, snr_eff_dB, p_draw, nvar_pilot);
            end

            clearvars pilot_real pilot_imag pilot_mask nvar_pilot ...
                      H_true h_tilde Ap delta_phi rainAtt_dB p_draw logp ...
                      snr_eff_dB snr_awgn pilotScale;

            sample_id = sample_id + 1;

        end % samples_per_esno loop
    end % esno_values loop

    fclose(fid);
    fprintf('\n  ✓ %s  →  %d rows written  →  %s\n\n', ...
        g.label, sample_id - 1, csvFile);

end % modcod_selection loop

fprintf('══════════════════════════════════════════════════════\n');
fprintf('  All MODCOD datasets complete.\n');
fprintf('  MODCODs     : %s\n', num2str(modcod_selection));
fprintf('  Es/No range : %.0f → %.0f dB  (step %g dB)\n', ...
    esno_values(1), esno_values(end), mean(diff(esno_values)));
fprintf('  Rows / file : %d  (%d levels × %d samples)\n', ...
    numel(esno_values)*samples_per_esno, numel(esno_values), samples_per_esno);
fprintf('  Output dir  : %s\n', outDir);
fprintf('  Ref pilots  : %s\n', refPilotFile);
fprintf('  Columns/row : %d meta + %d pilot = %d total  (per-MODCOD MAX_NP)\n', ...
    25, 3*max(cellfun(@(mc) geom(int32(mc)).Np, num2cell(modcod_selection))), ...
    25 + 3*max(cellfun(@(mc) geom(int32(mc)).Np, num2cell(modcod_selection))));
fprintf('══════════════════════════════════════════════════════\n');


%% 
%  HELPER FUNCTIONS  



function dfl = getDFL(modCod, fecFrame)
% GETDFL  Data field length (DFL) for a given MODCOD + FEC frame type.
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
    rxPilots  = rxFrame(pilotInd(:));
    totalPow  = mean(abs(rxPilots).^2);
    sigPow    = abs(mean(rxPilots .* conj(refPilots(:)))).^2;
    nVarEst   = max(0, totalPow - sigPow);
end