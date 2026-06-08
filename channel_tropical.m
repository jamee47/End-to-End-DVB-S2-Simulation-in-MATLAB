function chData = channel_tropical(txData, cfg, snr_dB, itu)

% Implements the full channel model from Awad et al. (2023), Eq. 8-11:
%
%   h_k = h_tilde * sqrt(b) * sqrt(b_max)           (Eq. 11)
%
% where:
%   h_tilde = A_p^(1/2) * exp(-j*phi)               (Eq. 8)
%   b       = G_s * beam_pattern^2                   (Eq. 10)
%   b_max   = (lam/4pi)^2 * G_R/(k_b*BW*T*d0^2)     (Eq. 9)
%
% ITU-R parameters (k, alpha, R_001, h_R) come from the itu struct
% produced by load_itu_params.m. Everything else is computed here.
%
% INPUT:
%   txData  — struct from tx_dvbs2.m
%   cfg     — configuration struct from main_generate_dataset.m
%   snr_dB  — SNR in dB
%   itu     — struct from load_itu_params.m
%
% OUTPUT:
%   chData.rxWaveform    — received waveform (h_tilde applied + AWGN)
%   chData.h_true        — full h_k (label for DL training / NMSE)
%   chData.h_tilde       — rain fading component (applied to waveform)
%   chData.b_beam        — satellite beam gain b
%   chData.b_max         — free space loss factor b_max
%   chData.rainAtten_dB  — instantaneous rain attenuation A_p (dB)
%   chData.A_001_dB      — worst-case rain attenuation (dB)
%   chData.snr_dB        — SNR used
%   chData.p_pct         — exceedance probability sampled (%)
% =========================================================================


required = {'k_V','alpha_V','k_H','alpha_H','R_001','h_R_km'};
for i = 1:length(required)
    if ~isfield(itu, required{i})
        error('[CH] Missing itu field: %s. Check load_itu_params.m', required{i});
    end
end


if isfield(cfg,'Polarization') && strcmpi(cfg.Polarization,'H')
    k_coef     = itu.k_H;
    alpha_coef = itu.alpha_H;
else
    k_coef     = itu.k_V;
    alpha_coef = itu.alpha_V;
end

% gamma_R = k * R_001^alpha   (dB/km)
gamma_R = k_coef * (itu.R_001 ^ alpha_coef);
fprintf('    [CH] gamma_R = %.4f dB/km\n', gamma_R);


% Full path calculation following ITU-R P.618-13
h_R       = itu.h_R_km;                 % rain height (km)
h_s       = cfg.Altitude_m / 1000;      % station altitude (km)
theta_deg = cfg.ElevationAngle;          % elevation angle (deg)
f_GHz     = cfg.Frequency_GHz;
lat       = itu.meta.lat;

% Slant path length through rain (km)
if theta_deg >= 5
    L_S = (h_R - h_s) / sind(theta_deg);
else
    L_S = 2*(h_R - h_s) / (sqrt(sind(theta_deg)^2 + ...
           2*(h_R - h_s)/8500) + sind(theta_deg));
end
L_G = L_S * cosd(theta_deg);            % horizontal projection (km)

% Horizontal reduction factor r_001 (P.618 Eq. 4)
r_001 = 1 / (1 + 0.78*sqrt(L_G*gamma_R/f_GHz) ...
              - 0.38*(1 - exp(-2*L_G)));

% Vertical adjustment factor v_001 (P.618 Eq. 5-7)
zeta = atand((h_R - h_s) / (L_G * r_001));
if zeta > theta_deg
    L_E = L_G * r_001 / cosd(theta_deg);
else
    L_E = (h_R - h_s) / sind(theta_deg);
end

chi = 36 - lat;
if abs(lat) < 36
    v_001 = 1 / (1 + sqrt(sind(theta_deg)) * ...
            (31*(1 - exp(-theta_deg/(1+chi))) * ...
             sqrt(L_E*gamma_R)/f_GHz^2 - 0.45));
else
    v_001 = 1 / (1 + sqrt(sind(theta_deg)) * ...
            (31*(1 - exp(-theta_deg)) * ...
             sqrt(L_E*gamma_R)/f_GHz^2 - 0.45));
end

A_001 = gamma_R * L_E * v_001;
fprintf('    [CH] A_001 = %.3f dB  (0.01%% exceedance)\n', A_001);


% A_p = A_001 * (p/0.01)^C
% C   = -(0.655 + 0.033*ln(p) - 0.045*ln(A_001) - 0.053*(1-p)*sin(el))
% p sampled from [0.1%, 10%] — realistic operating range

rng('shuffle');
p_pct  = 0.1 + (10 - 0.1) * rand(1);
C_exp  = -(0.655 + 0.033*log(p_pct) ...
             - 0.045*log(max(A_001, 0.01)) ...
             - 0.053*(1 - p_pct)*sind(theta_deg));
A_p_dB = A_001 * (p_pct / 0.01)^C_exp;
A_p_dB = max(A_p_dB,  0.1);    % floor: clear-sky residual
A_p_dB = min(A_p_dB, 20.0);    % cap: link unusable beyond this

fprintf('    [CH] A_p   = %.3f dB  (p = %.3f%% of year)\n', A_p_dB, p_pct);

f_Hz   = f_GHz * 1e9;
lambda = 3e8 / f_Hz;
d0     = 35788e3;                        % GEO altitude (m)
G_R    = 10^(cfg.G_r_dBi / 10);          % receiver gain (linear)
k_b    = 1.38e-23;                       % Boltzmann constant

b_max  = (lambda/(4*pi))^2 * (1/d0^2) * (G_R/(k_b * cfg.NoiseBW_Hz * cfg.NoiseTemp_K));

% b = G_s * [J1(u)/2u + 36*J3(u)/u^3]^2
% u = 2.07123 * sin(theta_offset) / sin(theta_3dB)

G_s   = 10^(cfg.G_s_dBi / 10);

if cfg.BeamOffset_deg < 1e-6
    beam_pattern_sq = 1.0;               % at beam centre: pattern = 1
else
    u = 2.07123 * sind(cfg.BeamOffset_deg) / sind(cfg.Theta3dB_deg);
    J1u = besselj(1, u);
    J3u = besselj(3, u);
    beam_pattern_sq = (J1u/(2*u) + 36*J3u/(u^3))^2;
end
b_beam = G_s * beam_pattern_sq;

fprintf('    [CH] b_max = %.4e | b_beam = %.4e\n', b_max, b_beam);

% h_tilde = A_p^(1/2) * exp(-j*phi),   phi ~ Uniform(0, 2*pi)
% Time-varying: resampled every frame call.

A_p_lin = 10^(-A_p_dB / 20);
phi     = 2 * pi * rand(1);
h_tilde = A_p_lin * exp(-1j * phi);

fprintf('    [CH] |h_tilde| = %.6f  angle = %.4f rad\n', ...
    abs(h_tilde), angle(h_tilde));


% Stored as h_true — this is the DL estimator training label.
% Only h_tilde is applied to the waveform (b and b_max are in link budget).

h_k = h_tilde * sqrt(b_beam) * sqrt(b_max);
fprintf('    [CH] |h_k|     = %.4e\n', abs(h_k));


fadedWaveform = txData.waveform * h_tilde;
rxWaveform    = awgn(fadedWaveform, snr_dB, 'measured');

% Pack output
chData.rxWaveform    = rxWaveform;
chData.h_true        = h_k;
chData.h_tilde       = h_tilde;
chData.b_beam        = b_beam;
chData.b_max         = b_max;
chData.rainAtten_dB  = A_p_dB;
chData.A_001_dB      = A_001;
chData.snr_dB        = snr_dB;
chData.p_pct         = p_pct;


end