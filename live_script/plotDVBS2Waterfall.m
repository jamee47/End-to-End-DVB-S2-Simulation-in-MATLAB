function plotDVBS2Waterfall(esno, berResults, perResults, varargin)
%PLOTDVBS2WATERFALL  Plot smoothed DVB-S2 BER/PER waterfall curves.
%
%   plotDVBS2Waterfall(esno, berResults, perResults)
%   plotDVBS2Waterfall(esno, berResults, perResults, 'Title', 'My Title')
%   plotDVBS2Waterfall(esno, berResults, perResults, 'Window', 5)
%   plotDVBS2Waterfall(esno, berResults, perResults, 'SaveFig', true)
%
%   The function applies three post-processing stages to convert noisy
%   Monte-Carlo BER/PER estimates into a clean waterfall curve:
%
%   Stage 1 — Validity filter
%       Remove NaN, zero, and saturated (= 1.0) points. Zero BER points
%       are artifacts of insufficient bit errors, not true zero BER.
%       Saturated points (= 1.0) indicate all frames lost.
%
%   Stage 2 — Log-domain moving average
%       Smooth in log10 space (not linear) so that noise near the
%       waterfall cliff — where BER spans several decades — is treated
%       proportionally. A symmetric window of ±Window/2 points is used.
%
%   Stage 3 — Monotone decreasing enforcement (isotonic regression)
%       A true BER/PER waterfall is strictly non-increasing with Es/No.
%       Any upward jags after smoothing are statistical artifacts and are
%       corrected by a cumulative-minimum scan from right to left (high
%       Es/No → low Es/No), which preserves the cliff location while
%       removing spurious bumps.

% ── Parse optional name-value arguments ───────────────────────────────────
p = inputParser;
addParameter(p, 'Title',    '', @ischar);
addParameter(p, 'Window',   5,  @(x) isnumeric(x) && x >= 1);
addParameter(p, 'SaveFig',  false, @islogical);
addParameter(p, 'FileName', 'dvbs2_waterfall', @ischar);
parse(p, varargin{:});

titleStr  = p.Results.Title;
winLen    = p.Results.Window;
doSave    = p.Results.SaveFig;
fileName  = p.Results.FileName;

% ── Input validation ───────────────────────────────────────────────────────
esno        = esno(:)';
berResults  = berResults(:)';
perResults  = perResults(:)';

assert(length(esno) == length(berResults), 'esno and berResults must be same length.');
assert(length(esno) == length(perResults), 'esno and perResults must be same length.');

% ── Process BER ────────────────────────────────────────────────────────────
[esno_ber, ber_raw, ber_smooth] = processWaterfall(esno, berResults, winLen);

% ── Process PER ────────────────────────────────────────────────────────────
[esno_per, per_raw, per_smooth] = processWaterfall(esno, perResults, winLen);

% ── Plot ───────────────────────────────────────────────────────────────────
figure('Color','k','Position',[100 100 860 540]);
ax = axes('Color','k','XColor','w','YColor','w','GridColor',[0.3 0.3 0.3], ...
    'MinorGridColor',[0.2 0.2 0.2],'GridAlpha',0.6,'MinorGridAlpha',0.4);
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
set(ax,'YScale','log','FontSize',11,'LineWidth',0.8);

% Raw data (faded markers only — shows where the noise came from)
if ~isempty(esno_ber)
    semilogy(ax, esno_ber, ber_raw, 'o', ...
        'Color',[0.27 0.51 0.96 0.35], ...   % faded blue
        'MarkerFaceColor','none', ...
        'MarkerSize', 5, ...
        'HandleVisibility','off');
end
if ~isempty(esno_per)
    semilogy(ax, esno_per, per_raw, 's', ...
        'Color',[0.96 0.30 0.30 0.35], ...   % faded red
        'MarkerFaceColor','none', ...
        'MarkerSize', 5, ...
        'HandleVisibility','off');
end

% Smoothed curves (solid line + filled markers)
if ~isempty(esno_ber)
    semilogy(ax, esno_ber, ber_smooth, '-o', ...
        'Color'   , [0.27 0.51 0.96], ...
        'MarkerFaceColor', [0.27 0.51 0.96], ...
        'MarkerSize', 5, ...
        'LineWidth', 2.0, ...
        'DisplayName', 'BER (smoothed)');
end
if ~isempty(esno_per)
    semilogy(ax, esno_per, per_smooth, '-s', ...
        'Color'   , [0.96 0.30 0.30], ...
        'MarkerFaceColor', [0.96 0.30 0.30], ...
        'MarkerSize', 5, ...
        'LineWidth', 2.0, ...
        'DisplayName', 'PER (smoothed)');
end

% Reference line: QEF threshold (BER = 10^-7 per DVB-S2 standard)
yline(ax, 1e-7, '--', 'QEF  (10^{-7})', ...
    'Color', [0.8 0.8 0.2], ...
    'LineWidth', 1.2, ...
    'LabelHorizontalAlignment','left', ...
    'HandleVisibility','off');

% Axes labels and title
xlabel(ax, 'E_s/N_0  (dB)', 'Color','w', 'FontSize', 12);
ylabel(ax, 'Error Rate',     'Color','w', 'FontSize', 12);

if isempty(titleStr)
    titleStr = 'DVB-S2 BER & PER Waterfall';
end
title(ax, titleStr, 'Color','w', 'FontSize', 13, 'FontWeight','bold');

% Dynamic y-axis limits
allSmooth = [ber_smooth, per_smooth];
if ~isempty(allSmooth)
    yMin = max(1e-8, 10^(floor(log10(min(allSmooth))) - 1));
    yMax = 1.5;
    ylim(ax, [yMin yMax]);
end
xlim(ax, [min(esno)-0.1, max(esno)+0.1]);

% Legend
leg = legend(ax, 'show', 'Location', 'southwest', 'TextColor', 'w', ...
    'Color', [0.12 0.12 0.12], 'EdgeColor', [0.4 0.4 0.4], 'FontSize', 10);

% Annotation: smoothing parameters
annotStr = sprintf('Window: %d pts  |  Smoothing: log-domain MA + monotone', winLen);
annotation('textbox',[0.13 0.01 0.75 0.04], ...
    'String', annotStr, ...
    'Color', [0.6 0.6 0.6], ...
    'EdgeColor', 'none', ...
    'FontSize', 8, ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', 'none');

% Save
if doSave
    exportgraphics(gcf, [fileName '.png'], 'Resolution', 200, 'BackgroundColor', 'k');
    fprintf('Figure saved → %s.png\n', fileName);
end

end % main function


% ══════════════════════════════════════════════════════════════════════════
%  Internal helper: apply all three smoothing stages
% ══════════════════════════════════════════════════════════════════════════
function [esno_valid, raw_valid, smoothed] = processWaterfall(esno, data, winLen)
%PROCESSWATERFALL  Filter → log-MA smooth → monotone enforce.

    % ── Stage 1: validity filter ──────────────────────────────────────────
    % Remove NaN, zero (no errors seen), and 1.0 (all frames lost).
    valid = ~isnan(data) & data > 0 & data < 1.0;
    esno_valid = esno(valid);
    raw_valid  = data(valid);

    if isempty(raw_valid)
        smoothed = [];
        return;
    end

    % ── Stage 2: log-domain moving average ────────────────────────────────
    % Operate in log10 space: noise near the waterfall cliff spans many
    % decades and is better averaged proportionally.
    logData = log10(raw_valid);
    n       = length(logData);
    logSmooth = zeros(1, n);
    half    = floor(winLen / 2);

    for k = 1:n
        lo = max(1, k - half);
        hi = min(n, k + half);
        logSmooth(k) = mean(logData(lo:hi));
    end

    % ── Stage 3: monotone decreasing enforcement ──────────────────────────
    % BER/PER must be non-increasing as Es/No increases. Scan from the
    % highest Es/No point leftward and propagate the minimum seen so far.
    % This removes upward jags without shifting the waterfall cliff.
    logMono = logSmooth;
    for k = n-1 : -1 : 1
        % logMono(k) must be >= logMono(k+1)  (higher BER at lower Es/No)
        logMono(k) = max(logMono(k), logMono(k+1));
    end

    smoothed = 10.^logMono;

end