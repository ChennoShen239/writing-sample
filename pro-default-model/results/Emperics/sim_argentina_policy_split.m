%% Argentina policy split simulation using real GDP (2006Q1–2013Q4)
% - Before 2007Q1: both paths use baseline policy
% - From 2007Q1 onward: Path A = baseline; Path B = theta=100
% - Outputs: debt/GDP (annualized) series and spread series

clear; clc;

%% Paths and utilities
rootDir = fileparts(fileparts(mfilename('fullpath'))); % .../results
resultsDir = fullfile(rootDir);
addpath(resultsDir); % to use formatFigure.m, loadBinary.m

% Model directories
baseDir = fullfile(resultsDir, 'baseline');
thetaDir = fullfile(resultsDir, 'thetad100');

%% Load model objects (grids, policies, prices, params)
baseline = loadModel(baseDir);
thetad100 = loadModel(thetaDir);

%% Load Argentina GDP (expects Emperics/data1.xlsx with columns: date, rgdp)
dataPath = fullfile(resultsDir, 'Emperics', 'data1.xlsx');
T = readtable(dataPath, 'VariableNamingRule', 'preserve');

% Parse quarter strings like "1987Q1 [1987Q1]" to datetime
s = string(T.date);
key = extractBefore(s, ' ');
yr  = str2double(extractBefore(key, 'Q'));
qtr = str2double(extractAfter(key,  'Q'));
mo  = (qtr - 1) * 3 + 1;
t   = datetime(yr, mo, 1);

% Numeric GDP series (log levels)
y_raw = T.rgdp;
if ~isnumeric(y_raw); y_raw = str2double(string(y_raw)); end

% Sample window and baseline split date
sampleMask = (t >= datetime(2005,12,1)) & (t <= datetime(2009,1,31)) & ~isnan(y_raw);
t_series = t(sampleMask);
y_series = y_raw(sampleMask);
logy_series = log(y_series);

splitDate = datetime(2007,1,1);
splitIx = find(t_series >= splitDate, 1, 'first');

%% Map realized log GDP to model y-states (quantile-standardization)
% Standardize model yGrid in logs
logy_grid = log(baseline.yGrid(:));
mu_model = mean(logy_grid);
sd_model = std(logy_grid);
z_model = (logy_grid - mu_model) / sd_model;

% Standardize data over the entire window for stable mapping
mu_data = mean(logy_series);
sd_data = std(logy_series);
z_data = (logy_series - mu_data) / sd_data;

% Nearest z-score match to pick the y-index each quarter
y_idx = zeros(length(z_data), 1);
for i = 1:length(z_data)
    [~, y_idx(i)] = min(abs(z_model - z_data(i)));
end

% Model-implied quarterly output level for each quarter (model units)
y_model_series = baseline.yGrid(y_idx);

%% Initialize debt at 2006Q1 to approx 25% of annual GDP (cap to grid)
% b_over_annual_gdp = b / (y * 4). Target 0.25 => b0_target = 0.25 * y * 4
b0_target = 0.055 * y_model_series(1) * 4;
b0_target = min(max(b0_target, baseline.bGrid(1)), baseline.bGrid(end));
[~, b_idx0] = min(abs(baseline.bGrid - b0_target));

%% Simulate two paths: A=baseline, B=thetad100 (after split)
Tn = length(t_series);

b_idx_A = zeros(Tn,1); b_idx_B = zeros(Tn,1);
b_A = zeros(Tn,1);    b_B = zeros(Tn,1);
q_A = zeros(Tn,1);    q_B = zeros(Tn,1);
sp_A = zeros(Tn,1);   sp_B = zeros(Tn,1);

b_idx_A(1) = b_idx0; b_idx_B(1) = b_idx0;
b_A(1) = baseline.bGrid(b_idx0);
b_B(1) = baseline.bGrid(b_idx0);

for tIx = 1:Tn
    yix = y_idx(tIx);
    % Select model per branch and time
    mA = baseline;                        % Path A always baseline
    mB = baseline; if tIx >= splitIx, mB = thetad100; end % Path B switches after split

    % Current debt indices
    bixA = b_idx_A(tIx);
    bixB = b_idx_B(tIx);

    % Transition over b' using policy distribution; compute expected b' and expected q
    pA = squeeze(mA.bPol(yix, bixA, :)); pA = pA / max(sum(pA), eps);
    pB = squeeze(mB.bPol(yix, bixB, :)); pB = pB / max(sum(pB), eps);

    % Expected next debt level and price (use dot products to avoid implicit expansion)
    EbA = sum(pA .* mA.bGrid);
    EbB = sum(pB .* mB.bGrid);
    EqA = (pA') * (mA.q(yix, :)');
    EqB = (pB') * (mB.q(yix, :)');

    % Spreads (quarterly to annualized) using pricing identity
    sp_q_A = (mA.rf + mA.delta) ./ max(EqA, eps) - mA.delta - mA.rf;
    sp_q_B = (mB.rf + mB.delta) ./ max(EqB, eps) - mB.delta - mB.rf;
    sp_A(tIx) = 4 * sp_q_A;  % annualized
    sp_B(tIx) = 4 * sp_q_B;

    % Store current q as expected q
    q_A(tIx) = EqA;
    q_B(tIx) = EqB;

    % Advance debt index for next period using nearest grid to expected b'
    if tIx < Tn
        [~, b_idx_A(tIx+1)] = min(abs(mA.bGrid - EbA));
        [~, b_idx_B(tIx+1)] = min(abs(mB.bGrid - EbB));
        b_A(tIx+1) = mA.bGrid(b_idx_A(tIx+1));
        b_B(tIx+1) = mB.bGrid(b_idx_B(tIx+1));
    end
end

% Debt-to-annual-GDP ratios (in %)
debt_ratio_A = 100 * (b_A ./ (y_model_series * 4));
debt_ratio_B = 100 * (b_B ./ (y_model_series * 4));

%% Monte Carlo bands for 95% confidence intervals (sampling b' from policy)
Nmc = 2000; rng(123);
[Bmc_A, SPmc_A] = simulateBranchMC(baseline, baseline, y_idx, splitIx, b_idx0, Nmc);
[Bmc_B, SPmc_B] = simulateBranchMC(baseline, thetad100, y_idx, splitIx, b_idx0, Nmc);

% Debt ratios for MC paths (Nmc x Tn)
den = (y_model_series(:)' * 4); % 1 x Tn
debt_ratio_mc_A = 100 * bsxfun(@rdivide, Bmc_A, den);
debt_ratio_mc_B = 100 * bsxfun(@rdivide, Bmc_B, den);

% Quantiles across simulations (rows = [2.5 50 97.5], cols = time)
drA_q = prctile(debt_ratio_mc_A, [2.5 50 97.5], 1);
drB_q = prctile(debt_ratio_mc_B, [2.5 50 97.5], 1);
spA_q = prctile(100*SPmc_A, [2.5 50 97.5], 1);
spB_q = prctile(100*SPmc_B, [2.5 50 97.5], 1);

%% Plot: Debt ratio series
figure(1); clf; set(gcf, 'Position', [100,100,900,600]);
hold on; xn = datenum(t_series);
% Colors consistent with repo style
colors = [188 228 183; 64 171 92; 0 78 45] / 255;
% 95% CI shading (baseline and theta=100)
shadeCI(t_series, drA_q(1,:), drA_q(3,:), colors(2,:), 0.15);
shadeCI(t_series, drB_q(1,:), drB_q(3,:), colors(3,:), 0.15);
hDebtA = plot(xn, debt_ratio_A, '--', 'Color', colors(2,:), 'LineWidth', 2.5, 'DisplayName', 'Baseline debt/GDP');
hDebtB = plot(xn, debt_ratio_B, '-',  'Color', colors(3,:), 'LineWidth', 2.5, 'DisplayName', 'High $\theta$ debt/GDP');
xlim([datenum(datetime(2005,12,31)) datenum(datetime(2009,1,1))]);

xsplit = datenum(splitDate);

% 图1
xline(xsplit, ':k', '2007-Q1', 'Interpreter','latex','FontSize',15,'LineWidth',2);


xlabel('Quarter', 'Interpreter','latex'); datetick('x','yyyy-QQ','keeplimits');
ylabel('Debt / Annual GDP (\%)', 'Interpreter','latex');
% title('Debt Ratio: Baseline vs High $\theta$ (Argentina GDP path)', 'Interpreter','latex');
grid on;
% xlim([datenum(datetime(2006,1,1)) datenum(datetime(2013,12,31))]);
formatLocal();
legend([hDebtA hDebtB], 'Location','southoutside', 'Orientation','horizontal', 'Interpreter','latex', 'Box','off');

% Save
out1 = fullfile(resultsDir, 'Emperics', 'argentina_debt_series.pdf');
drawnow;
print(gcf, out1, '-dpdf', '-r300', '-painters');

%% Plot: Spread series (annualized, %)
figure(2); clf; set(gcf, 'Position', [100,460,900,600]);
hold on;
% 95% CI shading for spreads
shadeCI(t_series, spA_q(1,:), spA_q(3,:), colors(2,:), 0.15);
shadeCI(t_series, spB_q(1,:), spB_q(3,:), colors(3,:), 0.15);
hSpA = plot(xn, 100*sp_A, '--', 'Color', colors(2,:), 'LineWidth', 2.5, 'DisplayName', 'Baseline spreads');
hSpB = plot(xn, 100*sp_B, '-',  'Color', colors(3,:), 'LineWidth', 2.5, 'DisplayName', 'High $\theta$ spreads');

xline(xsplit, ':k', '2007-Q1', 'Interpreter','latex','FontSize', 15,'LineWidth',2);


xlim([datenum(datetime(2005,12,31)) datenum(datetime(2009,1,1))]);
xlabel('Quarter', 'Interpreter','latex'); datetick('x','yyyy-QQ','keeplimits');
ylabel('Spread (annualized, \%)', 'Interpreter','latex');
% title('Spreads: Baseline vs High $\theta$ (Argentina GDP path)', 'Interpreter','latex');
grid on;
% xlim([datenum(datetime(2006,1,1)) datenum(datetime(2013,12,31))]);
formatLocal();
legend([hSpA hSpB], 'Location','southoutside', 'Orientation','horizontal', 'Interpreter','latex', 'Box','off');

out2 = fullfile(resultsDir, 'Emperics', 'argentina_spread_series.pdf');
drawnow;
print(gcf, out2, '-dpdf', '-r300', '-painters');

fprintf('\nSaved:\n  %s\n  %s\n', out1, out2);

%% Local loader for model objects
function M = loadModel(dirPath)
    % Read params
    p = dlmread(fullfile(dirPath, 'par.dat'));
    ySz = p(1); bSz = p(2); crra = p(3); rf = p(4); delta = p(5); %#ok<NASGU>

    % Grids
    M.yGrid = loadBinary(fullfile(dirPath, 'yGrid.bin'), 'float64', [ySz, 1]);
    M.bGrid = loadBinary(fullfile(dirPath, 'bGrid.bin'), 'float64', [bSz, 1]);

    % Objects
    M.q    = loadBinary(fullfile(dirPath, 'q.bin'),    'float64', [ySz, bSz]);
    M.dPol = loadBinary(fullfile(dirPath, 'dPol.bin'), 'float64', [ySz, bSz]);
    M.bPol = loadBinary(fullfile(dirPath, 'bPol.bin'), 'float64', [ySz, bSz, bSz]);
    M.EbPr = loadBinary(fullfile(dirPath, 'EbPr.bin'), 'float64', [ySz, bSz]);
    M.ySz = ySz; M.bSz = bSz; M.rf = rf; M.delta = delta; M.crra = crra;
end

%% Monte Carlo simulator for branch policies (returns Nmc x T arrays)
function [Bmc, SPmc] = simulateBranchMC(Mpre, Mpost, y_idx, splitIx, b_idx0, Nmc)
    Tn = length(y_idx);
    Bidx = zeros(Nmc, Tn);
    Bmc  = zeros(Nmc, Tn);
    SPmc = nan(Nmc, Tn);
    Bidx(:,1) = b_idx0;
    Bmc(:,1)  = Mpre.bGrid(b_idx0);
    min_q = 1e-6;
    for tIx = 1:Tn
        yix = y_idx(tIx);
        M = Mpre; if tIx >= splitIx, M = Mpost; end
        for n = 1:Nmc
            bix = Bidx(n,tIx);
            p = squeeze(M.bPol(yix, bix, :));
            s = sum(p);
            if s > 0
                p = p / s;
                u = rand;
                c = cumsum(p);
                bpr = find(c >= u, 1, 'first');
                if isempty(bpr); bpr = bix; end
            else
                % Fallback to nearest expected b' if degenerate distribution
                bpr_val = M.EbPr(yix, bix);
                [~, bpr] = min(abs(M.bGrid - bpr_val));
            end
            q_t = M.q(yix, bpr);
            sp_q = (M.rf + M.delta) / max(q_t, min_q) - M.delta - M.rf;
            SPmc(n,tIx) = 4 * sp_q;
            if tIx < Tn
                Bidx(n,tIx+1) = bpr;
                Bmc(n,tIx+1)  = M.bGrid(bpr);
            end
        end
    end
end

%% Shaded confidence band helper
function shadeCI(tvec, lb, ub, rgb, alpha)
    % Convert to numeric time for robust fill
    x = datenum(tvec(:));   % T x 1
    lb = lb(:);             % T x 1
    ub = ub(:);             % T x 1
    px = [x; flipud(x)];    % 2T x 1
    py = [lb; flipud(ub)];  % 2T x 1
    h = fill(px, py, rgb, 'EdgeColor','none'); %#ok<NASGU>
    set(h, 'FaceAlpha', alpha, 'HandleVisibility','off');
    datetick('x','yyyy-QQ','keeplimits');
end

%% Local formatting (self-contained)
function formatLocal()
    ax = gca;
    set(ax, 'FontSize', 15, 'TickLabelInterpreter', 'latex');
    set(findall(gcf,'Type','Line'), 'LineWidth', 2.5);
    set(gcf,'Color','w');
end
