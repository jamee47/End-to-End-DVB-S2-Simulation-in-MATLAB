%% Generating DVB-S2 waveform
% DVB-S2 configuration
wavegen = dvbs2WaveformGenerator;
wavegen.StreamFormat = "TS";
wavegen.NumInputStreams = 1;
wavegen.FECFrame = "normal";
wavegen.MODCOD = 1;
wavegen.DFL = 15928;
wavegen.HasPilots = 0;
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

Fs = 1e+06; 								 % Specify the sample rate of the waveform in Hz

% Filtering:

%% Visualize
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


