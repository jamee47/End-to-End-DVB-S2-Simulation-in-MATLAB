%% Generating DVB-S2 waveform
% DVB-S2 configuration
wavegen = dvbs2WaveformGenerator;
wavegen.StreamFormat = "TS";
wavegen.NumInputStreams = 1;
wavegen.FECFrame = "normal";
wavegen.MODCOD = 1;
wavegen.DFL = 15928;
wavegen.HasPilots = 1; % <--- CHANGED TO 1 TO ENABLE PILOTS
wavegen.RolloffFactor = 0.35;
wavegen.FilterSpanInSymbols = 10;
wavegen.SamplesPerSymbol = 4;

% input bit source:
NumPLFrames = 1;
pn = comm.PNSequence('Polynomial', 'x9+x5+1', 'InitialConditions', [zeros(1, 8) 1],'VariableSizeOutput',true,'MaximumOutputSize',[64800*NumPLFrames 1]);

if strcmp(wavegen.StreamFormat,"TS")
    syncBits = [0;1;0;0;0;1;1;1];
    pktLen = 1496;
    data =  cell(1, wavegen.NumInputStreams);
    for i = 1:wavegen.NumInputStreams
        numPkts = wavegen.MinNumPackets(i)*NumPLFrames;
        reset(pn);
        rawpkts = pn(pktLen * numPkts);
        txRawPkts = reshape(rawpkts, pktLen, numPkts);
        txPkts = [repmat(syncBits, 1, numPkts); txRawPkts];
        data{i} = txPkts(:);
    end
elseif strcmp(wavegen.StreamFormat,"GS") && any(wavegen.UPL>0)
    syncBits = [0;0;0;0;0;0;0;0];
    data =  cell(1, wavegen.NumInputStreams);
    for i = 1:wavegen.NumInputStreams
        numPkts = wavegen.MinNumPackets(i)*NumPLFrames;
        reset(pn);
        if length(wavegen.UPL) == 1
            pktLen = wavegen.UPL-8;
        else
            pktLen = wavegen.UPL(i)-8;
        end
        rawpkts = pn(pktLen * numPkts);
        txRawPkts = reshape(rawpkts, pktLen, numPkts);
        txPkts = [repmat(syncBits, 1, numPkts); txRawPkts];
        data{i} = txPkts(:);
    end
else
    data =  cell(1, wavegen.NumInputStreams);
    for i = 1:wavegen.NumInputStreams
        reset(pn);
        if length(wavegen.DFL) == 1
            rawpkts = pn(wavegen.DFL* NumPLFrames);
        else
            rawpkts = pn(wavegen.DFL(i)* NumPLFrames);
        end
        data{i} = rawpkts;
    end
end
in = data;

% Generation
waveform = [wavegen(in);flushFilter(wavegen)];

%% --- NEW CODE: Extract Symbols and Visualize Pilots ---

% 1. Matched Filter & Downsample to get back to symbol rate
txFilter = comm.RaisedCosineReceiveFilter( ...
    'RolloffFactor', wavegen.RolloffFactor, ...
    'FilterSpanInSymbols', wavegen.FilterSpanInSymbols, ...
    'InputSamplesPerSymbol', wavegen.SamplesPerSymbol, ...
    'DecimationFactor', wavegen.SamplesPerSymbol);

txSym_raw = txFilter(waveform);
filterDelay = wavegen.FilterSpanInSymbols / 2;
txSymbols = txSym_raw(filterDelay+1:end); % Remove filter delay transient
release(txFilter);

% 2. Calculate Pilot Indices (DVB-S2 Standard Structure)
% PLHEADER = 90 symbols. Pilots are 36 symbols inserted every 16 slots (1440 symbols).
PLHeaderLen = 90;
slotSize = 90;
slotsPerPilot = 16;
pilotBlockLen = 36;

symbolsPerDataChunk = slotsPerPilot * slotSize; % 1440 symbols
pilotSpacing = pilotBlockLen + symbolsPerDataChunk; % 1476 symbols

totalSymbols = length(txSymbols);
pilotIndices = [];

% Start looking for pilots after the first PLHEADER and first data chunk
pos = PLHeaderLen + symbolsPerDataChunk + 1; 

while pos + pilotBlockLen - 1 <= totalSymbols
    % Append the 36 indices for this pilot block
    pilotIndices = [pilotIndices, pos : (pos + pilotBlockLen - 1)]; 
    % Jump to the next pilot block
    pos = pos + pilotSpacing; 
end

% 3. Visualize the Real Part of the Symbols
figure('Name', 'DVB-S2 Frame Structure - Pilot Visualization', 'Position', [100, 100, 900, 400]);
% Plot all symbols in light gray
plot(1:totalSymbols, real(txSymbols), '.', 'Color', [0.7 0.7 0.7]);
hold on;
% Overlay pilot symbols in red
plot(pilotIndices, real(txSymbols(pilotIndices)), 'r.', 'MarkerSize', 10);

title('Transmitted DVB-S2 Symbols (Real Part)');
xlabel('Symbol Index');
ylabel('Amplitude (Real)');
legend('Data / Header Symbols', 'Pilot Blocks');
grid on;

% Zoom in on the first few pilot blocks to see them clearly
xlim([1, 6000]); 
ylim([-1.5 1.5]);

%% Original Visualizations
Fs = 1e+06; % Specify the sample rate of the waveform in Hz

% Spectrum Analyzer
spectrum = spectrumAnalyzer(SampleRate=Fs, AveragingMethod='exponential', ForgettingFactor=1);
spectrum(waveform);
release(spectrum);

% Constellation Diagram
constel = comm.ConstellationDiagram('ColorFading', true, ...
    'ShowTrajectory', 0, ...
    'SamplesPerSymbol', 4, ...
    'ShowReferenceConstellation', false);    

constel(waveform);
release(constel);

isPilot = false(size(txSymbols));
isPilot(pilotIndices) = true;

dataSymbols = txSymbols(~isPilot);
pilotSymbols = txSymbols(isPilot);

% 2. Create the Constellation Figure
figure('Name', 'DVB-S2 Constellation with Pilots', 'Position', [150, 150, 600, 600]);

% Plot Data Symbols in Yellow/Gold (similar to your screenshot)
scatter(real(dataSymbols), imag(dataSymbols), 10, [0.9290, 0.6940, 0.1250], 'filled', 'MarkerFaceAlpha', 0.3);
hold on;

% Plot Pilot Symbols in Red (slightly larger so they pop out)
scatter(real(pilotSymbols), imag(pilotSymbols), 30, 'r', 'filled');

% Formatting the plot to look like a standard constellation scope
title('DVB-S2 Constellation Diagram (QPSK)');
xlabel('In-phase Amplitude');
ylabel('Quadrature Amplitude');
legend('Data Symbols', 'Pilot Symbols', 'Location', 'best');
grid on;
axis square; % Keeps the I and Q axes proportionally equal

% Dynamically set limits based on the signal amplitude
maxAmp = max(abs([real(txSymbols); imag(txSymbols)])) * 1.2;
xlim([-maxAmp maxAmp]);
ylim([-maxAmp maxAmp]);

%% --- NEW: Extract Pilot Blocks and Export to CSV ---

% 1. Extract all pilot symbols using the indices we found earlier
allPilotSymbols = txSymbols(pilotIndices);

% 2. Reshape into a matrix: [numPilotBlocks x 36]
% MATLAB's reshape works column-wise, so we reshape to [36 x numBlocks] and transpose (.')
numBlocks = length(pilotIndices) / pilotBlockLen;
pilotMatrix = reshape(allPilotSymbols, pilotBlockLen, numBlocks).'; 

% 3. Split into Real and Imaginary parts for ML-friendly CSV export
% This creates a matrix of size [numBlocks x 72]
pilotMatrixReal = real(pilotMatrix);
pilotMatrixImag = imag(pilotMatrix);
exportMatrix = [pilotMatrixReal, pilotMatrixImag];

% 4. Create descriptive column headers for the CSV
headers = cell(1, pilotBlockLen * 2);
for i = 1:pilotBlockLen
    headers{i} = sprintf('Re_Sym_%02d', i);            % Columns 1-36: Real parts
    headers{i + pilotBlockLen} = sprintf('Im_Sym_%02d', i); % Columns 37-72: Imaginary parts
end

% 5. Convert the matrix to a Table and Export
pilotTable = array2table(exportMatrix, 'VariableNames', headers);
csvFilename = 'extracted_pilot_blocks.csv';
writetable(pilotTable, csvFilename);

% Display confirmation in the command window
fprintf('Successfully extracted %d pilot blocks.\n', numBlocks);
fprintf('Exported matrix of size [%d x %d] to %s\n', height(pilotTable), width(pilotTable), csvFilename);