% 定义颜色（保持原有配色方案）
colors = [188 228 183;     % 蓝色
    64 171 92; % 橙色
    0 78 45]; % 黄色

colors = colors / 255;

filename = 'data1';
try
    T = readtable(filename, 'VariableNamingRule', 'preserve');
catch
    % 旧版本 MATLAB 没有 VariableNamingRule 就回退
    T = readtable(filename);
end

% 输入：T（155x2 table），变量为 T.date（如 '1987Q1 [1987Q1]'）与 T.rgdp（double）

% 1) 解析季度时间
s = string(T.date);                         % e.g., "1987Q1 [1987Q1]"
key = extractBefore(s, ' ');                % -> "1987Q1"
yr  = str2double(extractBefore(key, 'Q'));  % -> 1987
q   = str2double(extractAfter(key,  'Q'));  % -> 1..4
mo  = (q - 1) * 3 + 1;                      % 1,4,7,10
t   = datetime(yr, mo, 1);                  % 季度起始月

% 2) 数值清洗（以防 rgdp 有字符或 '..'）
y = T.rgdp;
if ~isnumeric(y); y = str2double(string(y)); end

% 3) 选择 2006Q1–2013Q4 区间
mask = (t >= datetime(2006,1,1)) & (t <= datetime(2013,12,31)) & ~isnan(y);
t_sel = t(mask);
y_sel = y(mask);
y_sel = log(y_sel);
% 4) 绘图
figure;
plot(t_sel, y_sel, 'LineWidth', 1);
xlabel('Quarter');
ylabel('Log Real GDP (constant 2010 US$, millions, SA)');
title('Argentina: Log Quarterly Real GDP (2006Q1–2013Q4)');
grid on; set(gca, 'FontSize', 12);
formatFigure();

function formatFigure(hasLegendTitle, legendTitleText)
    % 格式化图像
    ax = gca; % Get current axes handle
    ax.FontSize = 18;
    ax.TickLabelInterpreter = 'latex';


    % Create the legend AFTER adjusting axes position
    % You might need to adjust FontSize for specific figures (e.g., Figures 2 & 4 had 10)
    lgd = legend('Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 3, ...
        'Interpreter', 'latex',  'Box', 'off');

    % Apply legend title if provided (for Figure 8)
    if nargin > 0 && hasLegendTitle
        title(lgd, legendTitleText);
    end

    % set(ax, 'Box', 'on', 'LineWidth', 1.5); % Uncomment if you want a box around the plot area itself
end
