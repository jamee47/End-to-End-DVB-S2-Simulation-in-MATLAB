function export_dataset_csv(allResults, cfg)
% =========================================================================
% EXPORT_DATASET_CSV — Save frame-level dataset to CSV
% =========================================================================
% Each row = one complete frame.
% Columns = [pb00_Re_y...pbN_Im_hLS | Re_hk | Im_hk |
%            Re_htilde | Im_htilde | MODCOD | SNR_dB |
%            RainAtten_dB | BER_LS | NMSE_LS]
%
% Reshape in Python:
%   X = df[feature_cols].values.reshape(-1, num_pilot_blocks, 6)
%   y = df[['Re_hk','Im_hk']].values
% =========================================================================

modNames      = fieldnames(allResults);
combinedTable = [];
headerBuilt   = false;
combinedHdr   = {};

for m = 1:length(modNames)
    name = modNames{m};
    res  = allResults.(name);

    datasetMatrix  = res.datasetMatrix;
    numPilotBlocks = res.numPilotBlocks;
    N              = size(datasetMatrix, 1);
    numFeatureCols = numPilotBlocks * 6;

    fprintf('    [CSV] %s: %d frames x %d pilot blocks x 6 features\n', ...
        name, N, numPilotBlocks);

    % Build column headers
    if ~headerBuilt
        hdr = {};
        for pb = 0:numPilotBlocks-1
            tag = sprintf('pb%02d', pb);
            hdr{end+1} = [tag '_Re_y'];
            hdr{end+1} = [tag '_Re_p'];
            hdr{end+1} = [tag '_Re_hLS'];
            hdr{end+1} = [tag '_Im_y'];
            hdr{end+1} = [tag '_Im_p'];
            hdr{end+1} = [tag '_Im_hLS'];
        end
        hdr{end+1} = 'Re_hk';
        hdr{end+1} = 'Im_hk';
        hdr{end+1} = 'Re_htilde';
        hdr{end+1} = 'Im_htilde';
        hdr{end+1} = 'SNR_dB';
        hdr{end+1} = 'RainAtten_dB';
        hdr{end+1} = 'BER_LS';
        hdr{end+1} = 'NMSE_LS';
        hdr{end+1} = 'BER_MMSE';
        hdr{end+1} = 'NMSE_MMSE';
        combinedHdr = hdr;
        headerBuilt = true;
    end

    T     = array2table(datasetMatrix, 'VariableNames', combinedHdr);
    fname = fullfile(cfg.OutputDir, sprintf('dataset_%s.csv', name));
    writetable(T, fname);
    fprintf('    [CSV] Saved: %s  [%d rows x %d cols]\n', fname, N, width(T));

    if isempty(combinedTable)
        combinedTable = datasetMatrix;
    else
        combinedTable = [combinedTable; datasetMatrix]; %#ok<AGROW>
    end
end

% Combined CSV
if ~isempty(combinedTable)
    T_all  = array2table(combinedTable, 'VariableNames', combinedHdr);
    fAll   = fullfile(cfg.OutputDir, 'dataset_combined.csv');
    writetable(T_all, fAll);
    fprintf('    [CSV] Combined: %s  [%d rows x %d cols]\n', ...
        fAll, height(T_all), width(T_all));
end

% Metadata JSON for Python
numPB = allResults.(modNames{1}).numPilotBlocks;
meta.num_pilot_blocks = numPB;
meta.num_features     = 6;
meta.label_cols       = {'Re_hk','Im_hk'};
meta.reshape_cmd      = sprintf('X.reshape(-1, %d, 6)', numPB);
meta.total_rows       = size(combinedTable, 1);
meta.feature_description = 'pb{N}_{Re|Im}_{y|p|hLS}';

metaFile = fullfile(cfg.OutputDir, 'dataset_metadata.json');
fid = fopen(metaFile, 'w');
fprintf(fid, '%s', jsonencode(meta));
fclose(fid);
fprintf('    [CSV] Metadata: %s\n', metaFile);
end
