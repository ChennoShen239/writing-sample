clc;
clear;

% Load your simulation results and grid data
% Make sure these files are in your current MATLAB path
load results.mat;    % Should contain bGrid, etc.
load simulation.mat; % Should contain bSimIx, spSim, cSim, gdpSim, tbSim, valid

% -------------------------------------------------------------------------
% 1. Prepare Data
% -------------------------------------------------------------------------

% Convert debt index to actual debt values
bSim = bGrid(bSimIx);
% bPrSim = bGrid(bPrSimIx); % bPrSim is usually not directly used for reporting standard moments

% Ensure 'valid' variable exists and is a logical array
% This is crucial for selecting the relevant simulation periods (e.g., excluding warm-up or default periods)
if ~exist('valid', 'var')
    % --- IMPORTANT: You need to define 'valid' based on your simulation logic ---
    % Example: If your Fortran simulation outputs from t=300 onwards,
    % and you want to exclude periods of default (dSimIx == 1).
    % Make sure dSimIx is loaded from simulation.mat or from your tab file.

    % For this example, let's assume valid periods are from 300 to the end, and not in default.
    % You will likely need to adjust 'warmUpPeriods' based on your Fortran output
    % and ensure 'dSimIx' (default indicator) is correctly loaded/extracted if you use it for validity.
    warmUpPeriods = 299; % Fortran code output loop starts from 300
    simLength = length(spSim);
    valid = false(1, simLength);
    valid(warmUpPeriods + 1 : end) = true;

    % If dSimIx is available and you want to exclude default periods from statistics:
    % Make sure dSimIx is loaded/created if not already in simulation.mat
    % if exist('dSimIx', 'var')
    %     valid = valid & (dSimIx == 0);
    % else
    %     disp('Warning: ''dSimIx'' variable not found. Cannot filter out default periods for validity.');
    % end
    disp('Warning: ''valid'' variable was not preloaded. A default ''valid'' array has been created.');
    disp('Please ensure its logic matches your desired analysis periods (e.g., excluding initial warm-up, default periods).');
end


% Logarithm of GDP and Consumption (common for calculating standard deviations and correlations)
loggdp = log(gdpSim(valid));
logc = log(cSim(valid));

% Trade Balance to GDP ratio
tby = tbSim(valid) ./ gdpSim(valid);

% -------------------------------------------------------------------------
% 2. Calculate and Print Annualized Statistics
% -------------------------------------------------------------------------

fprintf('--------------------------------------------\n');
fprintf('Model Moments (Annualized where appropriate)\n');
fprintf('--------------------------------------------\n');

% 1. Mean Debt to GDP Ratio (Annualized)
% bSim is end-of-quarter debt stock. gdpSim is quarterly GDP flow.
% To express debt as a ratio of *annual* GDP, we annualize quarterly GDP by multiplying by 4.
fprintf("Mean Debt to GDP (annualized)   %10.2f %%\n", 100.0 * mean(bSim(valid) ./ (gdpSim(valid) * 4)));

% 2. Mean Spread (Annualized)
% spSim represents a quarterly spread. To annualize, simply multiply by 4.
fprintf("Mean Spread (annualized)        %10.2f %%\n", 100.0 * mean(spSim(valid) * 4));

% 3. Standard Deviation of Spread (Annualized)
% Similarly, for spreads (rates), standard deviation is often linearly annualized by multiplying by 4.
fprintf("Std Spread (annualized)         %10.2f %%\n", 100.0 * std(spSim(valid) * 4));

% 4. Standard Deviation of Log Consumption (Annualized)
% For standard deviations of log variables (which approximate growth rates),
% we annualize by multiplying by sqrt(4) = 2, assuming quarterly growth rates are uncorrelated.
fprintf("Std log C (annualized)          %10.2f %%\n", 100.0 * std(logc) * 2);

% 5. Standard Deviation of Log GDP (Annualized)
% Annualized using the same logic as log consumption.
fprintf("Std log GDP (annualized)        %10.2f %%\n", 100.0 * std(loggdp) * 2);

% 6. Correlation: Spread and GDP (Quarterly)
% Correlation coefficients are dimensionless and are generally NOT annualized.
% We report the correlation between the quarterly series directly.
fprintf("Corr Sp, GDP (quarterly)        %10.2f %%\n", 100.0 * corr(spSim(valid), loggdp));

% 7. Correlation: Trade Balance/GDP and GDP (Quarterly)
% Correlation coefficients are not annualized.
fprintf("Corr TB/GDP, GDP (quarterly)    %10.2f %%\n", 100.0 * corr(tby, loggdp));

fprintf('--------------------------------------------\n');
