%% DVB-S2 Simulation in a Tropical Channel with ITU-R P.618 Rain Attenuation
clc; clearvars; close all;

%% Simulation Parameters
% Waveform parameters
esno = 30;                       % Es/N0 in dB
numFrames = 1;                   % Number of frames to simulate
FECframe = 'normal';             % FEC frame type
MODCOD = 1;                      % Modulation and Coding scheme index
RolloffFactor = 0.35;            % Rolloff factor for pulse shaping
FilterSpanInSymbols = 10;        % Filter span in symbols
SPS = 2;                         % Samples per symbol
HasPilots = true;                % Pilot insertion flag
DFL = getDFL(MODCOD, FECframe);  % Data Field Length based on parameters
filterDelay = FilterSpanInSymbols/2; % Filter delay

% Channel parameters: Tropical channel
chanParams = struct();
chanParams.Frequency_GHz = 20;           % Operating frequency in GHz
chanParams.Latitude = 5.17;              % Latitude in degrees
chanParams.Altitude_m = 57;              % Altitude in meters
chanParams.ElevationAngle = 45;          % Elevation angle in degrees
chanParams.Polarization = 'V';           % Polarization: 'V' or 'H'

% ITU rain parameters for tropical channel (Penang-like case)
chanParams.R001_mmph = 130;              % Rain rate R0.01 in mm/h
chanParams.IsothermHeight_km = 4.5;      % 0°C isotherm height in km
pRain = 0.01;                            % Percentage of time (rain attenuation exceedance)

% System parameters
simParams.G_s_dBi = 52.0;                % Satellite transmit gain in dBi
simParams.G_r_dBi = 41.7;                % Receive antenna gain in dBi
simParams.NoiseBW_Hz = 50e6;             % Noise bandwidth in Hz
simParams.NoiseTemp_K = 207;             % Noise temperature in Kelvin
simParams.Theta3dB_deg = 0.4;            % 3dB beamwidth in degrees
simParams.BeamOffset_deg = 0.0;          % Beam offset in degrees

%% Transmitter: Waveform Generation
% Create DVB-S2 waveform generator object
wavegen = dvbs2WaveformGenerator;
wavegen.StreamFormat = 'TS';
wavegen.NumInputStreams = 1;
wavegen.FECFrame = FECframe;
wavegen.MODCOD = MODCOD;
wavegen.HasPilots = HasPilots;
wavegen.RolloffFactor = RolloffFactor;
wavegen.FilterSpanInSymbols = FilterSpanInSymbols;
wavegen.SamplesPerSymbol = SPS;
wavegen.DFL = DFL;

% Generate input bits: sync bits + random PN sequence
syncBits = [0; 1; 0; 0; 0; 1; 1; 1];
pktLen = 1496;
pn = comm.PNSequence('Polynomial', 'x9+x5+1', ...
    'InitialConditions', [zeros(1,8) 1], ...
    'VariableSizeOutput', true, ...
    'MaximumOutputSize', [64800 * numFrames, 1]);
numPkts = wavegen.MinNumPackets(1) * numFrames;
reset(pn);
rawpkts = pn(pktLen * numPkts);
txRawPkts = reshape(rawpkts, pktLen, numPkts);
txPkts = [repmat(syncBits, 1, numPkts); txRawPkts];
inputBits = txPkts(:);
data = {inputBits};

% Generate waveform
txWaveform = [wavegen(data); flushFilter(wavegen)];

% Plot transmitted constellation
constelTx = comm.ConstellationDiagram( ...
    'Title', 'Transmitted Data', ...
    'SamplesPerSymbol', wavegen.SamplesPerSymbol, ...
    'ShowReferenceConstellation', false);
constelTx(txWaveform);
release(constelTx);

% Filter transmitted waveform
txFilt = comm.RaisedCosineReceiveFilter( ...
    'RolloffFactor', RolloffFactor, ...
    'FilterSpanInSymbols', FilterSpanInSymbols, ...
    'InputSamplesPerSymbol', SPS, ...
    'DecimationFactor', SPS);
txSym_raw = txFilt(txWaveform);
% FINE-TUNED: Apply symmetric truncation to match structural boundaries
txSymbols = txSym_raw(filterDelay + 1:end - filterDelay);
release(txFilt);

%% DVB-S2 Physical Layer Parameters
[modOrder, codeRate, cwLen] = satcom.internal.dvbs.getS2PHYParams(MODCOD, FECframe);
dataLen = cwLen / log2(modOrder);
slotLen = 90;
pilotBlkFreq = 16;
pilotSymLen = 36;

% Calculate number of pilot blocks
numPilotBlks = floor(dataLen / (slotLen * pilotBlkFreq));
if mod(dataLen, (slotLen*pilotBlkFreq)) == 0
    numPilotBlks = numPilotBlks - 1;
end
pilotLen = numPilotBlks * pilotSymLen;
plFrameSize = dataLen + pilotLen + slotLen;

% Display summary
fprintf(['\n           SUMMARY\n', ...
    'MODCOD                  : %d\n', ...
    'Modulation Order        : %d\n', ...
    'Code Rate               : %.5f\n', ...
    'LDPC Codeword Length    : %d bits\n', ...
    'Data Symbols/Frame      : %d\n', ...
    'Pilot Blocks/Frame      : %d (%d symbols each)\n', ...
    'PL Frame Size           : %d symbols\n'], ...
    MODCOD, modOrder, codeRate, cwLen, dataLen, ...
    numPilotBlks, pilotSymLen, plFrameSize);

%% Pilot Indices and Reference Pilots
% Generate scrambling sequence
plScrambIntSeq = satcom.internal.dvbs.plScramblingIntegerSequence(0);
cMap = [1; 1j; -1; -1j];
cSeq = cMap(plScrambIntSeq + 1);

% Pilot block indices
[~, pilotInd_raw] = satcom.internal.dvbs.pilotBlock(numPilotBlks);
pilotInd = pilotInd_raw + slotLen;

% Reference pilots
refPilots = (1 + 1j) / sqrt(2) .* cSeq(pilotInd_raw);
fprintf('Pilot indices range: [%d, %d]\n', min(pilotInd), max(pilotInd));
fprintf('refPilots first 3: %.3f%+.3fj  %.3f%+.3fj  %.3f%+.3fj\n', ...
    real(refPilots(1)), imag(refPilots(1)), ...
    real(refPilots(2)), imag(refPilots(2)), ...
    real(refPilots(3)), imag(refPilots(3)));

%% Channel: ITU-R P.618 Rain Attenuation
% Calculate rain attenuation Ap_dB
Ap_dB = calcRainAttenuationP618(chanParams, pRain);

% ==================== FIXED CHANNEL IMPLEMENTATION ====================
% 1. Find the power gain/loss factor Ap as defined in Eq. (7)
Ap = 10^(Ap_dB / 20); 
% 2. Generate random phase
phi = 2 * pi * rand(1); 
% 3. Calculate the fading coefficient h_tilde as defined in Eq. (8)
h_tilde = (Ap^(0.5)) * exp(-1i * phi); 
% 4. Apply attenuation on the waveform, matching the physical link direction.
h_attenuation = (1 / (Ap^(0.5))) * exp(-1i * phi); 


% Free-space link parameters
c_light = 3e8;                    % speed of light in m/s
kb = 1.38064852e-23;              % Boltzmann constant
freq_hz = chanParams.Frequency_GHz * 1e9;
lambda = c_light / freq_hz;
d0 = 35788e3;                     % GEO orbit distance in meters
d = 0;                            % beam center assumption
Gr_linear = 10^(simParams.G_r_dBi / 10);
b_max = (lambda / (4 * pi))^2 * (1 / (d0^2 + d^2)) * (Gr_linear / (kb * simParams.NoiseBW_Hz * simParams.NoiseTemp_K));

% Satellite beam gain
Gs_linear = 10^(simParams.G_s_dBi / 10);
theta = deg2rad(simParams.BeamOffset_deg);
theta3dB = deg2rad(simParams.Theta3dB_deg);
u_val = 2.07123 * sin(theta) / sin(theta3dB);
if abs(u_val) < 1e-8
    beamPattern = 1;
else
    beamPattern = besselj(1, u_val) / (2 * u_val) + 36 * besselj(3, u_val) / (u_val^3);
end
b_gain = Gs_linear * abs(beamPattern)^2;

% Final channel coefficient applying the attenuation drop
H_true = h_attenuation * sqrt(b_gain) * sqrt(b_max);

% Apply channel to waveform
attWaveform = H_true * txWaveform;

% Add AWGN noise
snr = esno - 10 * log10(SPS);
snr_linear = 10^(snr / 10);
txPower = mean(abs(txWaveform).^2);
noiseVar = txPower / snr_linear;
noise = sqrt(noiseVar/2) * (randn(size(attWaveform)) + 1i*randn(size(attWaveform)));
fadedWaveform = attWaveform + noise;

% Display channel summary
fprintf('\n========== CHANNEL MODEL SUMMARY ==========\n');
fprintf('Rain percentage p          : %.4f %%\n', pRain);
fprintf('Rainfall rate R0.01        : %.2f mm/h\n', chanParams.R001_mmph);
fprintf('Rain attenuation Ap_dB     : %.4f dB\n', Ap_dB);
fprintf('Rain fading h_tilde        : %.4e %+.4ei\n', real(h_tilde), imag(h_tilde));
fprintf('Free-space/link b_max      : %.4e\n', b_max);
fprintf('Satellite beam gain b_gain : %.4e\n', b_gain);
fprintf('Final channel H_true       : %.4e %+.4ei\n', real(H_true), imag(H_true));
fprintf('SNR used                   : %.2f dB\n', snr);
fprintf('===========================================\n');

%% Receiver: Demodulation and Processing
rxWaveform = fadedWaveform;

% Plot received constellation
constelRx = comm.ConstellationDiagram( ...
    'Title', 'Received Data', ...
    'SamplesPerSymbol', 1, ...
    'ShowReferenceConstellation', false);

% Fixed downsampling via DecimationFactor
rxFilter = comm.RaisedCosineReceiveFilter( ...
    'RolloffFactor', wavegen.RolloffFactor, ...
    'FilterSpanInSymbols', wavegen.FilterSpanInSymbols, ...
    'InputSamplesPerSymbol', wavegen.SamplesPerSymbol, ...
    'DecimationFactor', wavegen.SamplesPerSymbol);
rxSym_raw = rxFilter(rxWaveform);

% FINE-TUNED: Apply symmetric truncation to guarantee frame boundary synchronization
filterDelayMatched = wavegen.FilterSpanInSymbols/2;
rxSymbols = rxSym_raw(filterDelayMatched + 1:end - filterDelayMatched);
release(rxFilter);

% Plot constellation
constelRx(rxSymbols);
release(constelRx);

%% LS Channel Estimation, PL Header Recovery, Bit Recovery
berList = [];
nmseLSList = [];
HlsList = [];

for stIdx = 1:numFrames
    startIdx = (stIdx - 1) * plFrameSize + 1;
    endIdx = startIdx + plFrameSize - 1;
    
    if endIdx > length(rxSymbols)
        fprintf('\nFrame %d skipped: not enough received symbols.\n', stIdx);
        break;
    end
    
    rxFrame = rxSymbols(startIdx:endIdx);
    fprintf('\nFrame %d\n', stIdx);
    fprintf('Symbols in RX frame : %d\n', length(rxFrame));
    
    %% LS Estimator
    pilotBlock = rxFrame(pilotInd);
    LS_est = pilotBlock ./ refPilots;
    LS_matrix = reshape(LS_est, pilotSymLen, numPilotBlks);
    LS_mean_blk = mean(LS_matrix, 1);
    H_ls = mean(LS_mean_blk);
    HlsList = [HlsList; H_ls];
    
    % NMSE calculation
    nmseLS = abs(H_true - H_ls)^2 / abs(H_true)^2;
    nmseLSList = [nmseLSList; nmseLS];
    fprintf('H_ls                    : %.4e %+.4ei\n', real(H_ls), imag(H_ls));
    fprintf('LS NMSE                 : %.4e\n', nmseLS);
    
    %% Equalize using LS estimate
    rxFrameEqualized = rxFrame ./ H_ls;
    
    % Plot equalized constellation
    figure;
    scatterplot(rxFrameEqualized);
    title(['LS Estimation - Frame ' num2str(stIdx)]);
    
    %% PL Header Recovery
    rxPLHeader = rxFrameEqualized(1:90);
    try
        phyParams = dvbsPLHeaderRecover(rxPLHeader, "Mode", "DVB-S2/S2X regular");
        M = phyParams.ModulationOrder;
        if M == 0
            R = [];
        else
            R = eval(phyParams.LDPCCodeIdentifier);
        end
        fecFrame = phyParams.FECFrameLength;
        pilotStat = phyParams.HasPilots;
        
        if M ~= modOrder || R ~= codeRate || fecFrame ~= cwLen || ~pilotStat
            fprintf("PL header decoding failed\n");
            isHeaderRecovered = false;
        else
            fprintf("PL Header successfully decoded\n");
            disp(phyParams);
            isHeaderRecovered = true;
        end
    catch ME
        fprintf('PL header recovery error: %s\n', ME.message);
        isHeaderRecovered = false;
    end
    
    %% Bit Recovery
    if ~isHeaderRecovered
        fprintf('Skipping dvbs2BitRecover because PL header was not recovered correctly.\n');
        ber = 1.0;
        berList = [berList; ber];
        fprintf('Bit Error Rate (BER): 1.0000e+00\n');
        continue;
    end
    
    % Noise variance estimation
    nVar = DVBS2NoiseVarEstimate(rxFrameEqualized, pilotInd, refPilots, false);
    fprintf('Estimated noise variance : %.4e\n', nVar);
    
    % Bit recovery
    try
        % FINE-TUNED: Wrapped within a cell array to match multi-stream physical interface structures
        [decBitsTemp, isFrameLost] = dvbs2BitRecover({rxFrameEqualized}, nVar, false);
        recoveredBits = decBitsTemp;
    catch ME
        fprintf('dvbs2BitRecover failed: %s\n', ME.message);
        fprintf('Frame counted as lost.\n');
        ber = 1.0;
        berList = [berList; ber];
        fprintf('Bit Error Rate (BER): 1.0000e+00\n');
        continue;
    end
    
    % Calculate BER
    originalBits = vertcat(data{:});
    originalBits = originalBits(:);
    recoveredBits = vertcat(recoveredBits{:});
    recoveredBits = recoveredBits(:);
    
    if isempty(recoveredBits) || isempty(originalBits) || isFrameLost
        fprintf("WARNING: Frame lost or output 0 bits.\n");
        ber = 1.0;
        fprintf('Bit Error Rate (BER): 1.0000e+00\n');
    else
        minLen = min(length(originalBits), length(recoveredBits));
        origAligned = originalBits(1:minLen);
        recAligned = recoveredBits(1:minLen);
        [numErrors, ber] = biterr(origAligned, recAligned);
        totalBits = minLen;
        fprintf('Total Bits Compared : %d\n', totalBits);
        fprintf('Total Bit Errors    : %d\n', numErrors);
        fprintf('Bit Error Rate (BER): %.4e\n', ber);
    end
    berList = [berList; ber];
end

%% Final Summary
fprintf('\n========== FINAL LS RESULT SUMMARY ==========\n');
if ~isempty(berList)
    fprintf('Average BER             : %.4e\n', mean(berList));
else
    fprintf('Average BER             : Not available\n');
end
if ~isempty(nmseLSList)
    fprintf('Average LS NMSE         : %.4e\n', mean(nmseLSList));
else
    fprintf('Average LS NMSE         : Not available\n');
end
fprintf('=============================================\n');

%% Spectrum Analysis
BW = 36e6;
symbolRate = BW / (1 + wavegen.RolloffFactor);
Fs = symbolRate * wavegen.SamplesPerSymbol;
spectrum = spectrumAnalyzer( ...
    'SampleRate', Fs, ...
    'AveragingMethod', "exponential", ...
    'ForgettingFactor', 1, ...
    'ChannelNames', ["Transmitted Signal", "Received Signal"], ...
    'ShowLegend', true);
minLen = min(length(txWaveform), length(rxWaveform));
spectrum([txWaveform(1:minLen), rxWaveform(1:minLen)]);
release(spectrum);

%% Helper Functions
function dfl = getDFL(modCod, fecFrame)
    % Get data field length based on modulation and FEC frame type
    if strcmp(fecFrame, "normal")
        nDefVal = [16008 21408 25728 32208 38688 43040 48408 51648 53840 57472 ...
            58192 38688 43040 48408 53840 57472 58192 43040 48408 51648 53840 ...
            57472 58192 48408 51648 53840 57472 58192] - 80;
    else
        nDefVal = [3072 5232 6312 7032 9552 10632 11712 12432 13152 14232 0 ...
            9552 10632 11712 13152 14232 0 10632 11712 12432 13152 14232 0 11712 ...
            12432 13152 14232 0] - 80;
    end
    dfl = nDefVal(modCod);
end

function nVarEst = DVBS2NoiseVarEstimate(rxData, pilotInd, refPilots, normFlag)
    % Estimate noise variance
    if normFlag
        rxData = rxData / sqrt(mean(abs(rxData).^2));
    end
    rxPilots = rxData(pilotInd, :);
    pSigN = mean(abs(rxPilots).^2);
    pSig = abs(mean(rxPilots .* conj(refPilots)))^2;
    nVarEst = abs(mean(pSigN - pSig));
end

function Ap_dB = calcRainAttenuationP618(chanParams, p)
    % Calculate rain attenuation Ap_dB using ITU-R P.618
    f_GHz = chanParams.Frequency_GHz;
    theta_deg = chanParams.ElevationAngle;
    lat_deg = chanParams.Latitude;
    hs_km = chanParams.Altitude_m / 1000;
    R001 = chanParams.R001_mmph;
    h0_km = chanParams.IsothermHeight_km;
    hR_km = h0_km + 0.36;
    Re_km = 8500;
    
    if R001 <= 0 || hR_km <= hs_km
        Ap_dB = 0;
        return;
    end
    
    % Slant-path length
    if theta_deg >= 5
        Ls = (hR_km - hs_km) / sind(theta_deg);
    else
        Ls = 2 * (hR_km - hs_km) / (sqrt(sind(theta_deg)^2 + 2 * (hR_km - hs_km) / Re_km) + sind(theta_deg));
    end
    LG = Ls * cosd(theta_deg);
    
    % Specific attenuation coefficients
    [k, alpha] = ituP838_k_alpha(f_GHz, theta_deg, chanParams.Polarization);
    gamma_R = k * (R001^alpha);
    
    % Horizontal reduction factor
    r001 = 1 / (1 + 0.78 * sqrt((LG * gamma_R) / f_GHz) - 0.38 * (1 - exp(-2 * LG)));
    
    % Zeta for vertical adjustment
    zeta = atand((hR_km - hs_km) / (LG * r001));
    if zeta > theta_deg
        LR = (LG * r001) / cosd(theta_deg);
    else
        LR = (hR_km - hs_km) / sind(theta_deg);
    end
    
    % Additional correction based on latitude
    if abs(lat_deg) < 36
        chi = 36 - abs(lat_deg);
    else
        chi = 0;
    end
    
    v001 = 1 / (1 + sqrt(sind(theta_deg)) * (31 * (1 - exp(-theta_deg / (1 + chi))) * (sqrt(LR * gamma_R) / (f_GHz^2)) - 0.45));
    
    % Effective path length
    LE = LR * v001;
    
    % Attenuation exceeded for 0.01%
    A001 = gamma_R * LE;
    if abs(p - 0.01) < 1e-12
        Ap_dB = A001;
        return;
    end
    
    % For other percentages
    if p >= 1 || abs(lat_deg) >= 36
        beta = 0;
    elseif theta_deg >= 25
        beta = -0.005 * (abs(lat_deg) - 36);
    else
        beta = -0.005 * (abs(lat_deg) - 36) + 1.8 - 4.25 * sind(theta_deg);
    end
    
    if A001 <= 0
        Ap_dB = 0;
        return;
    end
    
    exponent = -(0.655 + 0.033 * log(p) - 0.045 * log(A001) - beta * (1 - p) * sind(theta_deg));
    Ap_dB = A001 * (p / 0.01)^exponent;
end

function [k, alpha] = ituP838_k_alpha(f_GHz, theta_deg, polarization)
    % Compute specific attenuation coefficients
    x = log10(f_GHz);
    
    % Coefficients for horizontal and vertical polarizations
    akH = [-5.33980, -0.35351, -0.23789, -0.94158];
    bkH = [-0.10008, 1.26970, 0.86036, 0.64552];
    ckH = [1.13098, 0.45400, 0.15354, 0.16817];
    mkH = -0.18961;
    ckH0 = 0.71147;
    
    akV = [-3.80595, -3.44965, -0.39902, 0.50167];
    bkV = [0.56934, -0.22911, 0.73042, 1.07319];
    ckV = [0.81061, 0.51059, 0.11899, 0.27195];
    mkV = -0.16398;
    ckV0 = 0.63297;
    
    logkH = sum(akH .* exp(-((x - bkH) ./ ckH).^2)) + mkH * x + ckH0;
    logkV = sum(akV .* exp(-((x - bkV) ./ ckV).^2)) + mkV * x + ckV0;
    kH = 10^logkH;
    kV = 10^logkV;
    
    % Alpha coefficients
    aaH = [-0.14318, 0.29591, 0.32177, -5.37610, 16.1721];
    baH = [1.82442, 0.77564, 0.63773, -0.96230, -3.29980];
    caH = [-0.55187, 0.19822, 0.13164, 1.47828, 3.43990];
    maH = 0.67849;
    caH0 = -1.95537;
    
    aaV = [-0.07771, 0.56727, -0.20238, -48.2991, 48.5833];
    baV = [2.33840, 0.95545, 1.14520, 0.791669, 0.791459];
    caV = [-0.76284, 0.54039, 0.26809, 0.116226, 0.116479];
    maV = -0.053739;
    caV0 = 0.83433;
    
    alphaH = sum(aaH .* exp(-((x - baH) ./ caH).^2)) + maH * x + caH0;
    alphaV = sum(aaV .* exp(-((x - baV) ./ caV).^2)) + maV * x + caV0;
    
    % Determine polarization
    pol = upper(string(polarization));
    if pol == "H"
        tau_deg = 0;
    elseif pol == "V"
        tau_deg = 90;
    else
        tau_deg = 45;
    end
    
    % Compute k and alpha based on polarization
    k = (kH + kV + (kH - kV) * cosd(theta_deg)^2 * cosd(2 * tau_deg)) / 2;
    alpha = (kH * alphaH + kV * alphaV + (kH * alphaH - kV * alphaV) * cosd(theta_deg)^2 * cosd(2 * tau_deg)) / (2 * k);
end