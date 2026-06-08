function pilotPos = detect_pilot_positions(txData, cfg)
% =========================================================================
% DETECT_PILOT_POSITIONS — Find actual pilot block positions in ONE frame
% =========================================================================
% Rather than assuming fixed pilot positions from the standard formula,
% this function DETECTS them from the actual generated waveform by
% correlating with the known pilot value.
%
% DVB-S2 pilots are un-modulated: all symbols = (1+j)/sqrt(2).
% A block of 36 consecutive symbols matching this pattern = pilot block.
%
% Called ONCE per MODCOD before the frame loop. The detected positions
% are then reused for every frame (pilot positions are deterministic for
% a given MODCOD and frame configuration).
%
% INPUT:
%   txData  — struct from tx_dvbs2_frame.m
%   cfg     — configuration struct
%
% OUTPUT:
%   pilotPos — vector of start indices (1-indexed) of each pilot block
%              in the symbol-rate sequence
% =========================================================================

txSymbols    = txData.txSymbols;
pilotValue   = txData.pilotValue;
pilotBlockLen = cfg.PilotBlockLen;   % 36 symbols
totalSyms    = length(txSymbols);

% Normalise symbols for correlation (remove amplitude)
txNorm = txSymbols / abs(pilotValue);

% Known pilot pattern: 36 symbols all equal to pilotValue (normalised = 1+j)
pilotPattern = ones(pilotBlockLen, 1);   % normalised pilot = 1+j / |1+j|

% Slide a window of length 36 across the symbol stream.
% Compute mean magnitude of (symbol - pilotValue) in each window.
% A pilot block will have near-zero error.

threshold = 0.3;   % max acceptable mean deviation (empirical)
pilotPos  = [];
minSpacing = cfg.SlotsPerPilot * cfg.SlotSize + pilotBlockLen - 5;  % ~1471

pos = cfg.PLHeaderLen + 1;   % start searching after PLHEADER

while pos + pilotBlockLen - 1 <= totalSyms
    window = txNorm(pos : pos + pilotBlockLen - 1);

    % Mean distance from normalised pilot value (1+j)/sqrt(2) normalised = 1
    deviation = mean(abs(window - (pilotValue/abs(pilotValue))));

    if deviation < threshold
        pilotPos(end+1) = pos; %#ok<AGROW>
        % Jump past this pilot block + next data chunk to avoid re-detection
        pos = pos + pilotBlockLen + minSpacing;
    else
        pos = pos + 1;
    end
end

if isempty(pilotPos)
    warning('[PILOT] No pilot blocks detected via correlation. Falling back to formula.');
    % Formula-based fallback (DVB-S2 standard)
    symbolsPerDataChunk = cfg.SlotsPerPilot * cfg.SlotSize;
    pos = cfg.PLHeaderLen + symbolsPerDataChunk + 1;
    while pos + pilotBlockLen - 1 <= totalSyms
        pilotPos(end+1) = pos; %#ok<AGROW>
        pos = pos + pilotBlockLen + symbolsPerDataChunk;
    end
end

fprintf('    [PILOT] Detected %d pilot blocks (first at sym %d, last at sym %d)\n', ...
    length(pilotPos), pilotPos(1), pilotPos(end));
end
