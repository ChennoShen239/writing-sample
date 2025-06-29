%% process.m - Script to load data, calculate and save moments
% This script should be placed in a directory containing simulation results
% (par.dat, sim.dat, and .bin files). It will calculate simulation moments
% and save them to 'moments.mat' in the same directory.
function process()
% The script should be run from the directory it is in.
dataDir = pwd;

fprintf('Processing data in %s...\n', dataDir);

% Load simulation data
try
    s = loadSimulationForProcessing(dataDir);
    % Calculate moments
    moments = calculateMomentsForProcessing(s);
    % Save moments
    save(fullfile(dataDir, 'moments.mat'), 'moments');
    fprintf('Moments saved to moments.mat in %s.\n', dataDir);
catch ME
    fprintf('Error processing data in %s: %s\n', dataDir, ME.message);
    rethrow(ME);
end
end


%% Local Helper Functions

function s = loadSimulationForProcessing(dataDir)
% 加载模拟数据并进行预处理
s = struct();

% 加载数据
try
    % Read parameters from par.dat
    param_file = fullfile(dataDir, 'par.dat');
    params = readmatrix(param_file);
    s.ySz = params(1);
    s.bSz = params(2);
    s.crra = params(3);
    s.rf = params(4);
    s.delta = params(5);

    % Read simulation data from sim.dat
    sim_file = fullfile(dataDir, 'sim.dat');
    opts = detectImportOptions(sim_file);
    opts.VariableNames = {'ySimIx', 'bSimIx', 'bPrSimIx', 'dSimIx', 'spSim', 'cSim', 'gdpSim', 'tbSim'};
    s.sim = readtable(sim_file, opts);

    % Load binary grid and policy files
    s.yGrid = read_bin_for_processing(fullfile(dataDir, 'yGrid.bin'), [s.ySz, 1]);
    s.bGrid = read_bin_for_processing(fullfile(dataDir, 'bGrid.bin'), [s.bSz, 1]);
    s.q = read_bin_for_processing(fullfile(dataDir, 'q.bin'), [s.ySz, s.bSz]);
    % These are not needed for moments calculation, so we don't load them.
    % s.V = read_bin_for_processing(fullfile(dataDir, 'V.bin'), [s.ySz, s.bSz]);
    % s.dPol = read_bin_for_processing(fullfile(dataDir, 'dPol.bin'), [s.ySz, s.bSz]);
    % s.bPol = read_bin_for_processing(fullfile(dataDir, 'bPol.bin'), [s.ySz, s.bSz, s.bSz]);
    % s.EbPr = read_bin_for_processing(fullfile(dataDir, 'EbPr.bin'), [s.ySz, s.bSz]);

catch ME
    fprintf('Error loading data for processing from %s\n', dataDir);
    rethrow(ME);
end

% 预处理模拟数据
s.T = height(s.sim);

% 使用循环从Grid和SimIx中正确构造模拟向量
s.bSim = zeros(s.T, 1);
for t = 1:s.T
    s.bSim(t) = s.bGrid(s.sim.bSimIx(t));
end

s.ySim = s.sim.gdpSim; % gdpSim is already the value, not index
s.spSim = s.sim.spSim;
s.cSim = s.sim.cSim;
s.tbSim = s.sim.tbSim;
s.defSim = s.sim.dSimIx;
s.gdpSim = s.sim.gdpSim;

% 排除过渡期
burnin = round(s.T * 0.1);
s.valid = (burnin + 1):s.T;
% In default periods, spread is set to a large negative number, filter that out
valid_indices_mask = s.spSim(s.valid) >= 0;
s.valid = s.valid(valid_indices_mask);
end


function moments = calculateMomentsForProcessing(s)
% this is calculateMoments from main.m
% 计算模拟数据的关键moments
valid_initial = s.valid;

% 从'q'策略函数中重新计算利差，以精确匹配原始方法
q_sim = zeros(length(valid_initial), 1);
y_ix_sim = s.sim.ySimIx(valid_initial);
b_pr_ix_sim = s.sim.bPrSimIx(valid_initial);

for i = 1:length(valid_initial)
    q_sim(i) = s.q(y_ix_sim(i), b_pr_ix_sim(i));
end

sp_quarterly = (s.rf + s.delta) ./ q_sim - s.delta - s.rf;

% 创建一个统一且健壮的过滤器，确保所有变量的数据对齐
final_mask = (q_sim > 0) & (sp_quarterly >= 0);

% -- 使用最终过滤器过滤所有相关的模拟向量 --
sp_quarterly_final = sp_quarterly(final_mask);

gdp = s.gdpSim(valid_initial);
gdp_final = gdp(final_mask);

b_stock = s.bSim(valid_initial);
b_stock_final = b_stock(final_mask);

c = s.cSim(valid_initial);
c_final = c(final_mask);

tb = s.tbSim(valid_initial);
tb_final = tb(final_mask);

% -- 从过滤后的向量计算所有moments --
log_gdp_final = log(gdp_final);
tb_over_gdp_final = tb_final ./ gdp_final;
b_over_annual_gdp_final = b_stock_final ./ (gdp_final * 4);

% 债务/GDP (相对于年度GDP, %)
moments.mean_debt_gdp = mean(b_over_annual_gdp_final) * 100;
moments.std_debt_gdp = std(b_over_annual_gdp_final) * 100;

% 利差 (年度化, %)
sp_annualized = sp_quarterly_final * 4;
moments.mean_spread = mean(sp_annualized) * 100;
moments.std_spread = std(sp_annualized) * 100;

% 波动率 (%)
moments.std_log_c = std(log(c_final)) * 2 * 100;
moments.std_log_gdp = std(log_gdp_final) * 2 * 100;

% 贸易差额/GDP (%)
moments.mean_tb_gdp = mean(tb_over_gdp_final) * 100;
moments.std_tb_gdp = std(tb_over_gdp_final) * 100;

% 相关性 (无单位)
if length(sp_quarterly_final) > 1 && length(log_gdp_final) > 1
    moments.corr_sp_gdp = corr(sp_quarterly_final, log_gdp_final);
else
    moments.corr_sp_gdp = NaN;
end

if length(tb_over_gdp_final) > 1 && length(log_gdp_final) > 1
    moments.corr_tb_gdp = corr(tb_over_gdp_final, log_gdp_final);
else
    moments.corr_tb_gdp = NaN;
end

if length(b_over_annual_gdp_final) > 1 && length(log_gdp_final) > 1
    moments.corr_debt_gdp = corr(b_over_annual_gdp_final, log_gdp_final);
else
    moments.corr_debt_gdp = NaN;
end

% 违约率 (%) - 在所有燃烧期后的周期内计算
burnin = round(s.T * 0.1);
post_burnin_indices = (burnin + 1):s.T;
moments.default_rate = mean(s.sim.dSimIx(post_burnin_indices)) * 100;
end

function data = read_bin_for_processing(file_path, dims)
% this is read_bin from main.m
f = fopen(file_path, 'r');
if f < 0
    error(['Cannot open file: ' file_path]);
end
data = fread(f, inf, 'float64');
fclose(f);
if prod(dims) ~= length(data)
    warning('Dimension mismatch when reading binary file %s. Expected %d elements, got %d.', file_path, prod(dims), length(data));
    % Fill with zeros if file is smaller than expected
    if length(data) < prod(dims)
        data_tmp = zeros(prod(dims), 1);
        data_tmp(1:length(data)) = data;
        data = data_tmp;
    else % Truncate if larger
        data = data(1:prod(dims));
    end
end

% Reshape the 1D vector read from file into the correct N-D matrix dimensions.
% Both Fortran and Matlab use column-major order, so a direct reshape is correct.
data = reshape(data, dims);
end