%% DVB-S2 BER vs Es/No — Symbol-Level Baseband Simulation
%
% Generates a proper waterfall BER curve by simulating entirely at the
% symbol level.  This avoids the waveform-generation / frame-sync /
% dvbs2BitRecover failure chain that caused the flat BER in the original
% code.
%
% The three bugs in the original code that produced flat lines:
%   1. rxSymbols(1:plFrameSize) never advanced the read pointer, so every
%      loop iteration re-processed the same first frame.
%   2. dvbs2BitRecover returned empty cells on header-decode failures,
%      forcing ber=1.0 every frame and hiding the SNR dependence.
%   3. sqrt(Ap_linear) applied the attenuation twice in dB.
%
% This script replaces that entire chain with a direct encode→mod→
% AWGN→demod→decode loop that cannot fail mid-frame.
%
% SUPPORTED MODCOD VALUES (pass any integer 1-28 except 11,17,23,28 for short)
%   1  QPSK  1/4     9  QPSK  8/9    17  8PSK  9/10   25  32APSK 8/9
%   2  QPSK  1/3    10  QPSK  9/10   18  16APSK 2/3   26  32APSK 9/10
%   3  QPSK  2/5    11  (reserved)   19  16APSK 3/4   ... etc.
%   4  QPSK  1/2    12  8PSK  3/5    20  16APSK 4/5
%   5  QPSK  3/5    13  8PSK  2/3    21  16APSK 5/6
%   6  QPSK  2/3    14  8PSK  3/4    22  16APSK 8/9
%   7  QPSK  3/4    15  8PSK  4/5    23  (reserved)
%   8  QPSK  4/5    16  8PSK  5/6    24  32APSK 3/4

clc; clearvars; close all;

%% =========================================================================
%  USER PARAMETERS
% ==========================================================================
MODCOD      = 16;              % DVB-S2 MODCOD index (1-28)
FECframe    = 'normal';       % 'normal' (64800 bit) or 'short' (16200 bit)
EsNodB_range = 0 : 0.5 : 8;  % Es/No sweep (dB) — bracket the threshold
numFrames   = 30;             % frames per SNR point (more = smoother curve)
maxLDPCIter = 50;             % LDPC decoder iterations
matFileName = 'dvbs2_ber_results.mat';
% ==========================================================================

%% ---- Resolve MODCOD → physical layer parameters ------------------------
[modOrder, codeRate, cwLen] = satcom.internal.dvbs.getS2PHYParams(MODCOD, FECframe);
bitsPerSym = log2(modOrder);

% LDPC matrices
if strcmp(FECframe,'normal')
    H = dvbs2ldpc(codeRate);
else
    H = dvbs2ldpc(codeRate, 'short');
end
encCfg = ldpcEncoderConfig(H);
decCfg = ldpcDecoderConfig(H);

% k_ldpc: systematic bits returned by ldpcDecode (= cols - rows of H)
k_ldpc = size(H,2) - size(H,1);

% Bit interleaver order (DVB-S2 Sec 5.3.3)
% Required for 8PSK and higher to decorrelate LDPC parity-bit errors
if modOrder > 4
    intrlvOrder = getInterleaverOrder(cwLen, bitsPerSym);
else
    intrlvOrder = [];   % QPSK: no interleaver
end

fprintf('\n=================================================================\n');
fprintf(' DVB-S2 Symbol-Level BER Sweep\n');
fprintf('=================================================================\n');
fprintf(' MODCOD %2d | M=%-2d | Rate=%.5f | Frame=%s | cwLen=%d\n', ...
        MODCOD, modOrder, codeRate, FECframe, cwLen);
fprintf(' k_ldpc=%d | bitsPerSym=%d | frames/SNR=%d\n', ...
        k_ldpc, bitsPerSym, numFrames);
fprintf(' SNR range : %.1f to %.1f dB (%d points)\n', ...
        EsNodB_range(1), EsNodB_range(end), numel(EsNodB_range));
fprintf('=================================================================\n\n');

%% ---- SNR sweep ----------------------------------------------------------
numSNRpts  = numel(EsNodB_range);
BER_demod  = zeros(1, numSNRpts);
BER_LDPC   = zeros(1, numSNRpts);

for snrIdx = 1:numSNRpts
    EsNodB  = EsNodB_range(snrIdx);

    % Complex noise variance for unit-power symbols:
    %   Es/N0 = 10^(EsNodB/10)  →  N0 = 1/10^(EsNodB/10)
    %   noise n ~ CN(0, N0): real & imag each N(0, N0/2)
    N0     = 10^(-EsNodB/10);
    sigma  = sqrt(N0/2);        % std dev per real/imag component

    errDemod = 0;  totDemod = 0;
    errLDPC  = 0;  totLDPC  = 0;

    for frIdx = 1:numFrames

        %% -- (1) GENERATE RANDOM MESSAGE BITS ----------------------------
        txMsgBits = randi([0 1], k_ldpc, 1);

        %% -- (2) LDPC ENCODE  --------------------------------------------
        txCodeword = ldpcEncode(txMsgBits, encCfg);   % length = cwLen

        %% -- (3) BIT INTERLEAVE (8PSK and above) -------------------------
        if ~isempty(intrlvOrder)
            txInterleaved = txCodeword(intrlvOrder);
        else
            txInterleaved = txCodeword;
        end

        %% -- (4) MODULATE ------------------------------------------------
        txSymbols = modulateDVBS2(txInterleaved, modOrder, codeRate);

        %% -- (5) AWGN CHANNEL --------------------------------------------
        %  Direct noise injection — no awgn() wrapper to avoid SNR
        %  reference ambiguity.
        noise     = sigma * (randn(size(txSymbols)) + 1j*randn(size(txSymbols)));
        rxSymbols = txSymbols + noise;

        %% -- (6) SOFT DEMODULATE → LLR -----------------------------------
        llr = demodDVBS2soft(rxSymbols, modOrder, codeRate, N0);
        llr = llr(:);   % ensure column vector, length = cwLen

        %% -- (7) PRE-LDPC BER (hard decision vs TX codeword) -------------
        rxHard   = double(llr < 0);
        errDemod = errDemod + sum(xor(rxHard, txInterleaved));
        totDemod = totDemod + cwLen;

        %% -- (8) BIT DE-INTERLEAVE LLR (8PSK and above) -----------------
        if ~isempty(intrlvOrder)
            llrDeint            = zeros(cwLen,1);
            llrDeint(intrlvOrder) = llr;
            llr = llrDeint;
        end

        %% -- (9) LDPC DECODE ---------------------------------------------
        try
            [decBits, ~] = ldpcDecode(llr, decCfg, maxLDPCIter, ...
                                      'DecisionType', 'hard');
            % decBits has length k_ldpc (systematic bits only)
        catch
            decBits = ones(k_ldpc,1,'int8');  % count as all errors
        end

        %% -- (10) POST-LDPC BER (vs original message bits) ---------------
        errLDPC = errLDPC + sum(xor(logical(decBits(1:k_ldpc)), ...
                                    logical(txMsgBits)));
        totLDPC = totLDPC + k_ldpc;

    end  % frame loop

    BER_demod(snrIdx) = errDemod / max(totDemod, 1);
    BER_LDPC(snrIdx)  = errLDPC  / max(totLDPC,  1);

    fprintf('Es/No=%5.1f dB | BER_demod=%.3e | BER_LDPC=%.3e\n', ...
            EsNodB, BER_demod(snrIdx), BER_LDPC(snrIdx));

end  % SNR loop

%% ---- Save results -------------------------------------------------------
snrdB_sweep = EsNodB_range;
save(matFileName, 'snrdB_sweep','BER_demod','BER_LDPC', ...
     'MODCOD','FECframe','modOrder','codeRate','cwLen','k_ldpc','numFrames');
fprintf('\nResults saved → %s\n', matFileName);

%% ---- Plot ---------------------------------------------------------------
plotWaterfall(matFileName);


%% =========================================================================
%  LOCAL FUNCTIONS
% ==========================================================================

function sym = modulateDVBS2(bits, modOrder, codeRate)
%MODULATEDVBS2  Map bits to complex symbols using DVB-S2 constellations.
    switch modOrder
        case 4   % QPSK: pi/4 offset, Gray coding (DVB-S2 Sec 5.4.1)
            sym = pskmod(bits, 4, pi/4, 'gray', 'InputType','bit');
        case 8   % 8PSK: DVB-S2 Table 9 mapping
            sym = pskmod(bits, 8, 0, 'gray', 'InputType','bit');
        case 16  % 16APSK
            sym = dvbsapskmod(bits, 16, 's2', codeRate, ...
                              'InputType','bit','UnitAveragePower',true);
        case 32  % 32APSK
            sym = dvbsapskmod(bits, 32, 's2', codeRate, ...
                              'InputType','bit','UnitAveragePower',true);
        otherwise
            error('Unsupported modulation order: %d', modOrder);
    end
end

function llr = demodDVBS2soft(sym, modOrder, codeRate, N0)
%DEMODDVBS2SOFT  Compute soft LLRs — must use same params as modulateDVBS2.
    switch modOrder
        case 4
            llr = pskdemod(sym, 4, pi/4, 'gray', ...
                           'OutputType','approxllr','NoiseVariance',N0);
        case 8
            llr = pskdemod(sym, 8, 0, 'gray', ...
                           'OutputType','approxllr','NoiseVariance',N0);
        case 16
            llr = dvbsapskdemod(sym, 16, 's2', codeRate, ...
                                'OutputType','approxllr','NoiseVar',N0, ...
                                'UnitAveragePower',true);
        case 32
            llr = dvbsapskdemod(sym, 32, 's2', codeRate, ...
                                'OutputType','approxllr','NoiseVar',N0, ...
                                'UnitAveragePower',true);
        otherwise
            error('Unsupported modulation order: %d', modOrder);
    end
    llr = llr(:);
end

function order = getInterleaverOrder(cwLen, bitsPerSym)
%GETINTERLEAVERORDER  DVB-S2 Sec 5.3.3 column-wise bit interleaver.
%  Writes bits row-by-row into a (cwLen/bitsPerSym) x bitsPerSym matrix
%  then reads column-by-column.
    nRows  = cwLen / bitsPerSym;
    cols   = bitsPerSym;
    idx    = reshape(1:cwLen, cols, nRows)';  % fill row-by-row
    order  = idx(:);                          % read column-by-column
end

function plotWaterfall(matFile)
%PLOTWATERFALL  Load saved results and draw the SNR-vs-BER waterfall.
    S = load(matFile);

    b1 = S.BER_demod;  b1(b1 == 0) = NaN;
    b2 = S.BER_LDPC;   b2(b2 == 0) = NaN;

    figure('Name','DVB-S2 BER Waterfall','NumberTitle','off', ...
           'Position',[80 80 820 520]);

    semilogy(S.snrdB_sweep, b1, 'b--o','LineWidth',1.5,'MarkerSize',6, ...
             'DisplayName','Pre-LDPC (demod output)');
    hold on;
    semilogy(S.snrdB_sweep, b2, 'r-s','LineWidth',2,'MarkerSize',7, ...
             'DisplayName','Post-LDPC (end-to-end)');
    hold off;

    grid on; grid minor;
    xlabel('Es/No (dB)','FontSize',13);
    ylabel('Bit Error Rate (BER)','FontSize',13);
    title(sprintf('DVB-S2 BER Waterfall — MODCOD %d (M=%d, R=%.5f, %s frame)', ...
                  S.MODCOD, S.modOrder, S.codeRate, S.FECframe), ...
          'FontSize',12);
    legend('Location','southwest','FontSize',11);
    ylim([1e-6 1]);

    annotation('textbox',[0.13 0.01 0.6 0.05], ...
               'String', sprintf('Frames per SNR point: %d  |  max LDPC iter: %d', ...
               S.numFrames, 50), ...
               'EdgeColor','none','FontSize',8,'Color',[0.45 0.45 0.45]);
end