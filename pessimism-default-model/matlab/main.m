%% 主脚本：对比 baseline (thetaD=1) 和 high thetaD (thetaD=100) 的结果
clear; close all; clc;

% 获取当前脚本所在的目录
[scriptDir, ~, ~] = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end
cd(scriptDir);

%% 加载baseline数据
fprintf('Loading baseline data (thetaD=1)...\n');
try
    baselineDir = fullfile(scriptDir, 'baseline');
    baseline = loadSimulation(baselineDir);
    fprintf('Baseline data loaded successfully.\n');
catch ME
    fprintf('Failed to load baseline data: %s\n', ME.message);
    rethrow(ME);
end

%% 加载thetad100数据
fprintf('Loading high thetaD data (thetaD=100)...\n');
try
    thetad100Dir = fullfile(scriptDir, 'thetad100');
    thetad100 = loadSimulation(thetad100Dir);
    fprintf('High thetaD data loaded successfully.\n');
catch ME
    fprintf('Failed to load thetad100 data: %s\n', ME.message);
    rethrow(ME);
end

%% 加载theta10数据
fprintf('Loading medium thetaD data (thetaD=10)...\n');
try
    theta10Dir = fullfile(scriptDir, 'theta10');
    theta10 = loadSimulation(theta10Dir);
    fprintf('Medium thetaD data loaded successfully.\n');
catch ME
    fprintf('Failed to load theta10 data: %s\n', ME.message);
    rethrow(ME);
end

%% 计算并对比模拟moments
fprintf('\n=== SIMULATION MOMENTS COMPARISON ===\n');
original_dir = pwd;

try
    % 计算baseline的moments
    fprintf('\nCalculating baseline moments...\n');
    baseline_moments = calMoments(baseline);

    % 计算thetad100的moments
    fprintf('\nCalculating high thetaD moments...\n');
    thetad100_moments = calMoments(thetad100);

    % 计算theta10的moments
    fprintf('\nCalculating medium thetaD moments...\n');
    theta10_moments = calMoments(theta10);
catch ME
    cd(original_dir);
    rethrow(ME);
end

% 显示对比结果
fprintf('\n');
fprintf('%-25s %15s %15s %15s\n', 'Moment', 'Baseline', 'Med $\theta$', 'High $\theta$');
fprintf('%s\n', repmat('-', 1, 75));

moments_fields = {'mean_debt_gdp', 'std_debt_gdp', 'mean_spread', 'std_spread', ...
    'std_log_c', 'std_log_gdp', 'corr_sp_gdp', 'corr_tb_gdp', 'mean_tb_gdp', 'std_tb_gdp'};
moments_names = {'Mean Debt/GDP', 'Std Debt/GDP', 'Mean Spread', 'Std Spread', ...
    'Std log C', 'Std log GDP', 'Corr(Sp,GDP)', ...
    'Corr(TB/GDP,GDP)', 'Mean TB/GDP', 'Std TB/GDP'};

for i = 1:length(moments_fields)
    field = moments_fields{i};
    name = moments_names{i};

    % 检查字段是否存在
    if isfield(baseline_moments, field) && isfield(theta10_moments, field) && isfield(thetad100_moments, field)
        baseline_val = baseline_moments.(field);
        theta10_val = theta10_moments.(field);
        thetad100_val = thetad100_moments.(field);
        fprintf('%-25s %15.2f %15.2f %15.2f\n', name, 100*baseline_val, 100*theta10_val, 100*thetad100_val);
    else
        fprintf('%-25s %15s %15s %15s\n', name, 'N/A', 'N/A', 'N/A');
    end
end

%% 保存moments对比结果
fprintf('\nSaving moments comparison...\n');
moments_comparison = struct();
moments_comparison.baseline = baseline_moments;
moments_comparison.theta10 = theta10_moments;
moments_comparison.thetad100 = thetad100_moments;
save('moments_comparison.mat', 'moments_comparison');

plots;

%% Helper Functions

function s = loadSimulation(dataDir)
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
    s.yGrid = read_bin(fullfile(dataDir, 'yGrid.bin'), [s.ySz, 1]);
    s.bGrid = read_bin(fullfile(dataDir, 'bGrid.bin'), [s.bSz, 1]);
    s.V = read_bin(fullfile(dataDir, 'V.bin'), [s.ySz, s.bSz]);
    s.dPol = read_bin(fullfile(dataDir, 'dPol.bin'), [s.ySz, s.bSz]);
    s.q = read_bin(fullfile(dataDir, 'q.bin'), [s.ySz, s.bSz]);
    s.bPol = read_bin(fullfile(dataDir, 'bPol.bin'), [s.ySz, s.bSz, s.bSz]);
    s.EbPr = read_bin(fullfile(dataDir, 'EbPr.bin'), [s.ySz, s.bSz]);

catch ME
    fprintf('Error loading data from %s\n', dataDir);
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

K = 20;
N = 20;
s.valid = false([s.T, 1]);
for ix = K+N+1:s.T
    s.valid(ix) = (sum(s.defSim(ix-N:ix)) == 0);
end
end



function data = read_bin(file_path, dims)
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


function m = calMoments(s)
% ---------------------------------------------------------------------
% 这个函数现在是完全正确和规范的
% ---------------------------------------------------------------------

% 1. 准备数据
% 直接使用 loadSimulation 已经准备好的数据向量
valid_bSim   = s.bSim(s.valid);
valid_gdpSim = s.gdpSim(s.valid);
valid_spSim  = s.spSim(s.valid);
valid_cSim   = s.cSim(s.valid);
valid_tbSim  = s.tbSim(s.valid);

% 计算衍生变量
loggdp = log(valid_gdpSim);
logc   = log(valid_cSim);

% 计算每个季度的(年化)债务/GDP比率序列
% 使用括号让 (valid_gdpSim * 4) 更清晰地表示年化GDP
debt_gdp_ratio_annualized = valid_bSim ./ (valid_gdpSim * 4);

% 计算贸易差额与GDP比率序列
tby = valid_tbSim ./ valid_gdpSim;

% 2. 计算并存储各个moment
m = struct();

% 债务相关
m.mean_debt_gdp = mean(debt_gdp_ratio_annualized);
m.std_debt_gdp  = std(debt_gdp_ratio_annualized);

% 利差相关
m.mean_spread = mean(valid_spSim) * 4;
m.std_spread  = std(valid_spSim) * 4;

% 宏观波动相关
m.std_log_c = std(logc) * 2;
m.std_log_gdp = std(loggdp) * 2;

% 相关性
m.corr_sp_gdp = corr(valid_spSim, loggdp);
m.corr_tb_gdp = corr(tby, loggdp);

% 贸易平衡相关
m.mean_tb_gdp = mean(tby);
m.std_tb_gdp  = std(tby);

end







