function data = loadSimulation(dataDir)
%LOADSIMULATION 从指定目录加载模拟结果数据
%   data = loadSimulation(dataDir) 从 dataDir 目录加载所有结果文件
%   返回包含所有数据的结构体

if nargin < 1
    dataDir = './';
end

% 确保目录路径以 / 结尾
if ~endsWith(dataDir, '/')
    dataDir = [dataDir '/'];
end

fprintf('Attempting to load data from directory: %s\n', dataDir);
fullParamPath = fullfile(dataDir, 'parameters.tab');
fprintf('Attempting to read parameters from: %s\n', fullParamPath);
if ~exist(fullParamPath, 'file')
    error('Parameters file does not exist at the specified path: %s', fullParamPath);
end

fprintf('Loading data from: %s\n', dataDir);

try
    % 首先读取参数文件
    params = dlmread([dataDir 'parameters.tab']);
    ix = 1;
    ySz = params(ix); ix = ix + 1;
    bSz = params(ix); ix = ix + 1;
    crra = params(ix); ix = ix + 1;
    rf = params(ix); ix = ix + 1;
    delta = params(ix);

    % 使用loadBinary函数读取二进制文件
    data.yGrid = loadBinary([dataDir 'yGrid.bin'], 'float64', [ySz, 1]);
    data.yPi = loadBinary([dataDir 'yPi.bin'], 'float64', [ySz, ySz]);
    data.bGrid = loadBinary([dataDir 'bGrid.bin'], 'float64', [bSz, 1]);

    data.V = loadBinary([dataDir 'V.bin'], 'float64', [ySz, bSz]);
    data.Vr = loadBinary([dataDir 'Vr.bin'], 'float64', [ySz, bSz]);
    data.Vd = loadBinary([dataDir 'Vd.bin'], 'float64', [ySz, 1]);
    data.q = loadBinary([dataDir 'q.bin'], 'float64', [ySz, bSz]);
    data.dPol = loadBinary([dataDir 'dPol.bin'], 'float64', [ySz, bSz]);
    data.bPol = loadBinary([dataDir 'bPol.bin'], 'float64', [ySz, bSz, bSz]);
    data.EbPr = loadBinary([dataDir 'EbPr.bin'], 'float64', [ySz, bSz]);

    % 存储参数
    data.ySz = ySz;
    data.bSz = bSz;
    data.crra = crra;
    data.rf = rf;
    data.delta = delta;

    % 读取模拟数据（如果存在）
    simFile = [dataDir 'sim.tab'];
    if exist(simFile, 'file')
        simData = dlmread(simFile);
        ix = 1;
        data.ySimIx = simData(:, ix); ix = ix + 1;
        data.bSimIx = simData(:, ix); ix = ix + 1;
        data.bPrSimIx = simData(:, ix); ix = ix + 1;
        data.dSimIx = simData(:, ix); ix = ix + 1;
        data.spSim = simData(:, ix); ix = ix + 1;
        data.cSim = simData(:, ix); ix = ix + 1;
        data.gdpSim = simData(:, ix); ix = ix + 1;
        data.tbSim = simData(:, ix);

        % 转换spread
        data.spSim = (1 + data.spSim).^4 - 1;

        % 计算valid标记
        sz = size(data.spSim, 1);
        K = 20;
        N = 20;
        data.valid = false([sz, 1]);
        for ix = K+N+1:sz
            data.valid(ix) = (sum(data.dSimIx(ix-N:ix)) == 0);
        end

        fprintf('Simulation data loaded successfully.\n');
    else
        fprintf('No simulation data found.\n');
    end

    fprintf('Data loaded successfully from %s\n', dataDir);

catch ME
    fprintf('Error loading data from %s: %s\n', dataDir, ME.message);
    rethrow(ME);
end

end