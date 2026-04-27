function itu = load_itu_params(json_file)
% =========================================================================
% LOAD_ITU_PARAMS — Load ITU-R parameters from get_parameters.py JSON
% =========================================================================
% Reads the JSON output of get_parameters.py and returns a clean struct.
% This struct is passed directly into channel_tropical.m.
%
% Usage:
%   itu = load_itu_params('channel_params.json');
%   chData = channel_tropical(txData, cfg, snr_dB, itu);
%
% INPUT:
%   json_file — path to JSON file from get_parameters.py
%
% OUTPUT itu struct:
%   .k_V          — k coefficient, vertical pol.       [P.838]
%   .alpha_V      — alpha coefficient, vertical pol.   [P.838]
%   .k_H          — k coefficient, horizontal pol.     [P.838]
%   .alpha_H      — alpha coefficient, horizontal pol. [P.838]
%   .R_001        — rainfall rate at 0.01% exceed.(mm/h) [P.837]
%   .h_0_km       — 0 deg C isotherm height (km)       [P.839]
%   .h_R_km       — rain height h_R = h_0+0.36 (km)   [P.839]
%   .A_gas_dB     — gaseous attenuation (dB)           [P.676] optional
%   .A_cloud_dB   — cloud attenuation (dB)             [P.840] optional
%   .A_scint_dB   — scintillation attenuation (dB)     [P.618] optional
%   .meta         — location and link metadata sub-struct
% =========================================================================

%% Validate
if nargin < 1 || isempty(json_file)
    error('load_itu_params: json_file path is required.');
end
if ~isfile(json_file)
    error('load_itu_params: File not found: %s\nRun get_parameters.py first.', ...
        json_file);
end

%% Read and parse JSON
fprintf('[ITU] Loading: %s\n', json_file);
try
    mr = jsondecode(fileread(json_file));
    mr = mr.matlab_ready;
catch ME
    error('load_itu_params: JSON parse failed — %s', ME.message);
end

%% Extract fields
itu.k_V        = mr.k_V;
itu.alpha_V    = mr.alpha_V;
itu.k_H        = mr.k_H;
itu.alpha_H    = mr.alpha_H;
itu.R_001      = mr.R_001_mmh;
itu.h_0_km     = mr.h_0_km;
itu.h_R_km     = mr.h_R_km;
itu.A_gas_dB   = mr.A_gas_dB;      % 0.0 if not computed
itu.A_cloud_dB = mr.A_cloud_dB;    % 0.0 if not computed
itu.A_scint_dB = mr.A_scint_dB;    % 0.0 if not computed

%% Store metadata
meta           = jsondecode(fileread(json_file));
meta           = meta.metadata;
itu.meta.lat   = meta.latitude_deg;
itu.meta.lon   = meta.longitude_deg;
itu.meta.alt_m = meta.altitude_m;
itu.meta.freq  = meta.frequency_GHz;
itu.meta.el    = meta.elevation_deg;
itu.meta.pol   = meta.polarization;

%% Print summary
fprintf('[ITU] Location : (%.2fN, %.2fE)  Alt: %.0f m\n', ...
    itu.meta.lat, itu.meta.lon, itu.meta.alt_m);
fprintf('[ITU] Link     : %.1f GHz | El: %.1f deg | Pol: %s\n', ...
    itu.meta.freq, itu.meta.el, itu.meta.pol);
fprintf('[ITU] %-25s k=%.8f  alpha=%.8f\n', 'P.838 (V pol):', itu.k_V, itu.alpha_V);
fprintf('[ITU] %-25s k=%.8f  alpha=%.8f\n', 'P.838 (H pol):', itu.k_H, itu.alpha_H);
fprintf('[ITU] %-25s %.4f mm/h\n', 'R_001 (P.837):',  itu.R_001);
fprintf('[ITU] %-25s %.4f km\n',   'h_0   (P.839):',  itu.h_0_km);
fprintf('[ITU] %-25s %.4f km\n',   'h_R   (P.839):',  itu.h_R_km);
if mr.A_gas_dB > 0
    fprintf('[ITU] %-25s %.4f dB\n', 'A_gas  (P.676):', itu.A_gas_dB);
    fprintf('[ITU] %-25s %.4f dB\n', 'A_cloud(P.840):', itu.A_cloud_dB);
    fprintf('[ITU] %-25s %.4f dB\n', 'A_scint(P.618):', itu.A_scint_dB);
else
    fprintf('[ITU] Optional params (A_gas, A_cloud, A_scint): not in JSON.\n');
    fprintf('[ITU] Re-run get_parameters.py with --optional to include.\n');
end
fprintf('[ITU] Parameters loaded OK.\n\n');
end
