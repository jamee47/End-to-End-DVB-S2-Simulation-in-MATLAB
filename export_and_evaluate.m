function export_and_evaluate(allResults, cfg)

% EXPORT_DATASET_CSV — Save frame-level dataset to CSV for BLSTM/GRU

% Dataset structure (one row = one complete frame):
%
%   X: [NumFrames x (NumPilotBlocks * 6)] — flattened sequence
%      Columns: pb_00_Re_y, pb_00_Re_p, pb_00_Re_hLS,
%               pb_00_Im_y, pb_00_Im_p, pb_00_Im_hLS,
%               pb_01_Re_y, ... pb_N_Im_hLS
%
%   y: [NumFrames x 2] — channel coefficient label
%      Columns: label_Re_hk, label_Im_hk
%
%   Meta: frame_id, MODCOD, SNR_dB, RainAtten_dB
%
% In Python, reshape X back to (NumFrames, NumPilotBlocks, 6)
% before feeding to BLSTM/GRU:
%   X = df[feature_cols].values.reshape(N, num_pilot_blocks, 6)
%   y = df[['label_Re_hk','label_Im_hk']].values
%
% Why one row per frame:
%   BLSTM/GRU treat each row as one training sample.
%   The sequence of pilot blocks within a frame is the temporal
%   input — splitting them into separate rows destroys this structure.
% =========================================================================

modNames      = fieldnames(allResults);
combinedRows  = {};
combinedHdr   = {};
headerBuilt   = false;

for m = 1:length(modNames)
    name   = modNames{m};
    rxData = allResults.(name);

    if isempty(rxData.X_in) || rxData.NumPilotBlocks == 0
        fprintf('    [CSV] Skipping %s — no data.\n', name);
        continue;
    end

    % X_in is [NumPilotBlocks x 6] for the ENTIRE run (all frames merged)
    % We need to split it back into per-frame sequences.
    % NumPilotBlocks per frame is stored in rxData.PilotBlocksPerFrame
    % If not stored, estimate from total / NumFrames

    totalBlocks = rxData.NumPilotBlocks;

    if isfield(rxData, 'PilotBlocksPerFrame')
        pb_per_frame = rxData.PilotBlocksPerFrame;
    else
        pb_per_frame = floor(totalBlocks / cfg.NumFrames);
        warning('[CSV] PilotBlocksPerFrame not found, estimating as %d', ...
            pb_per_frame);
    end

    if pb_per_frame == 0
        fprintf('    [CSV] Skipping %s — pb_per_frame = 0.\n', name);
        continue;
    end

    % Number of complete frames we can build
    numCompleteFrames = floor(totalBlocks / pb_per_frame);
    fprintf('    [CSV] %s: %d frames x %d pilot blocks x 6 features\n', ...
        name, numCompleteFrames, pb_per_frame);

    % Build column headers once
    if ~headerBuilt
        hdr = {};
        for pb = 0:pb_per_frame-1
            tag = sprintf('pb%02d', pb);
            hdr{end+1} = [tag '_Re_y'];
            hdr{end+1} = [tag '_Re_p'];
            hdr{end+1} = [tag '_Re_hLS'];
            hdr{end+1} = [tag '_Im_y'];
            hdr{end+1} = [tag '_Im_p'];
            hdr{end+1} = [tag '_Im_hLS'];
        end
        hdr{end+1} = 'label_Re_hk';
        hdr{end+1} = 'label_Im_hk';
        hdr{end+1} = 'frame_id';
        hdr{end+1} = 'MODCOD';
        hdr{end+1} = 'SNR_dB';
        hdr{end+1} = 'RainAtten_dB';
        combinedHdr = hdr;
        headerBuilt = true;
    end

    % Build one row per frame
    frameRows = zeros(numCompleteFrames, pb_per_frame*6 + 4);

    for f = 1:numCompleteFrames
        % Extract pilot blocks belonging to this frame
        startIdx = (f-1)*pb_per_frame + 1;
        endIdx   =  f   *pb_per_frame;
        block    = rxData.X_in(startIdx:endIdx, :);  % [pb_per_frame x 6]

        % Flatten row-major: pb0_Re_y, pb0_Re_p, ..., pb0_Im_hLS, pb1_Re_y,...
        frameRows(f, 1:pb_per_frame*6) = reshape(block', 1, []);

        % Labels: Re(h_k) and Im(h_k) — same for whole frame
        frameRows(f, pb_per_frame*6+1) = real(rxData.h_true);
        frameRows(f, pb_per_frame*6+2) = imag(rxData.h_true);

        % Metadata
        frameRows(f, pb_per_frame*6+3) = f;
        frameRows(f, pb_per_frame*6+4) = rxData.MODCOD;
    end

    % Append SNR and rain columns
    SNR_col  = rxData.snr_dB      * ones(numCompleteFrames, 1);
    rain_col = rxData.rainAtten_dB * ones(numCompleteFrames, 1);
    frameRows = [frameRows, SNR_col, rain_col]; %#ok<AGROW>

    % Update header for the extra columns
    if length(combinedHdr) == pb_per_frame*6 + 4
        combinedHdr{end+1} = 'SNR_dB';
        combinedHdr{end+1} = 'RainAtten_dB';
    end

    % Save per-MODCOD CSV
    T     = array2table(frameRows, 'VariableNames', combinedHdr);
    fname = fullfile(cfg.OutputDir, sprintf('dataset_%s.csv', name));
    writetable(T, fname);
    fprintf('    [CSV] Saved: %s  [%d rows x %d cols]\n', ...
        fname, numCompleteFrames, width(T));

    combinedRows{end+1} = frameRows; %#ok<AGROW>
end

% Save combined CSV
if ~isempty(combinedRows)
    combined     = vertcat(combinedRows{:});
    T_combined   = array2table(combined, 'VariableNames', combinedHdr);
    combinedFile = fullfile(cfg.OutputDir, 'dataset_combined.csv');
    writetable(T_combined, combinedFile);
    fprintf('    [CSV] Combined: %s  [%d rows x %d cols]\n', ...
        combinedFile, height(T_combined), width(T_combined));

    % Save reshape instructions for Python
    meta.num_pilot_blocks = pb_per_frame;
    meta.num_features     = 6;
    meta.feature_names    = {'Re_y','Re_p','Re_hLS','Im_y','Im_p','Im_hLS'};
    meta.label_names      = {'label_Re_hk','label_Im_hk'};
    meta.reshape_note     = sprintf( ...
        'X = df[feature_cols].values.reshape(-1, %d, 6)', pb_per_frame);
    meta.total_rows       = height(T_combined);

    metaFile = fullfile(cfg.OutputDir, 'dataset_metadata.json');
    fid = fopen(metaFile, 'w');
    fprintf(fid, '%s', jsonencode(meta));
    fclose(fid);
    fprintf('    [CSV] Metadata: %s\n', metaFile);
end
end