%% 创建对比图
fprintf('\nCreating comparison plots...\n');

% 定义线条样式
lineWidth = 3;  % 更粗的线条
baselineStyle = '--';  % 虚线
theta10Style = '-.'; % 点划线
thetad100Style = '-';   % 实线

% 定义颜色（保持原有配色方案）
colors = [188 228 183;     % 蓝色
    64 171 92; % 橙色
    0 78 45]; % 黄色

colors = colors / 255;
% 选择不同的y状态索引（低、中、高产出）
fixYs = round(baseline.ySz/2) + [-10 0 10];

%% 图1：价值函数对比 (Value Function)
figure(1);
set(gcf, 'Position', [100, 100, 800, 600]);
hold on;
plot(baseline.bGrid, baseline.V(fixYs(1),:), baselineStyle, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline Low');
plot(baseline.bGrid, baseline.V(fixYs(2),:), baselineStyle, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline Med');
plot(baseline.bGrid, baseline.V(fixYs(3),:), baselineStyle, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline High');
plot(theta10.bGrid, theta10.V(fixYs(1),:), theta10Style, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ Low');
plot(theta10.bGrid, theta10.V(fixYs(2),:), theta10Style, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ Med');
plot(theta10.bGrid, theta10.V(fixYs(3),:), theta10Style, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ High');
plot(thetad100.bGrid, thetad100.V(fixYs(1),:), thetad100Style, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ Low');
plot(thetad100.bGrid, thetad100.V(fixYs(2),:), thetad100Style, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ Med');
plot(thetad100.bGrid, thetad100.V(fixYs(3),:), thetad100Style, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ High');
xlabel('$B$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$V(B,y)$', 'Interpreter', 'latex', 'FontSize', 12);
% IMPORTANT: Remove legend call from here. It will be handled by formatFigure().
% legend('Location', 'best', 'Interpreter', 'latex');
grid on;
xlim([0 0.7]); % 规定横轴为从0 到 0.7
formatFigure(); % Call formatFigure after all plotting is done

%% 图2：违约概率对比 (Default Probability)
figure(2);
set(gcf, 'Position', [150, 150, 800, 600]);
hold on;
plot(baseline.bGrid, baseline.dPol(fixYs(1),:), baselineStyle, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline Low');
plot(baseline.bGrid, baseline.dPol(fixYs(2),:), baselineStyle, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline Med');
plot(baseline.bGrid, baseline.dPol(fixYs(3),:), baselineStyle, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline High');
plot(theta10.bGrid, theta10.dPol(fixYs(1),:), theta10Style, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ Low');
plot(theta10.bGrid, theta10.dPol(fixYs(2),:), theta10Style, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ Med');
plot(theta10.bGrid, theta10.dPol(fixYs(3),:), theta10Style, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ High');
plot(thetad100.bGrid, thetad100.dPol(fixYs(1),:), thetad100Style, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ Low');
plot(thetad100.bGrid, thetad100.dPol(fixYs(2),:), thetad100Style, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ Med');
plot(thetad100.bGrid, thetad100.dPol(fixYs(3),:), thetad100Style, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ High');
xlabel('$B$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$P(\mathrm{default})$', 'Interpreter', 'latex', 'FontSize', 12);
% IMPORTANT: Remove legend call from here.
% legend('Location', 'best', 'Interpreter', 'latex');
grid on;
xlim([0 0.7]); % 规定横轴为从0 到 0.7
formatFigure();

%% 图3：债券价格对比 (Bond Price)
figure(3);
set(gcf, 'Position', [200, 200, 800, 600]);
hold on;
plot(baseline.bGrid, baseline.q(fixYs(1),:), baselineStyle, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline Low');
plot(baseline.bGrid, baseline.q(fixYs(2),:), baselineStyle, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline Med');
plot(baseline.bGrid, baseline.q(fixYs(3),:), baselineStyle, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline High');
plot(theta10.bGrid, theta10.q(fixYs(1),:), theta10Style, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ Low');
plot(theta10.bGrid, theta10.q(fixYs(2),:), theta10Style, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ Med');
plot(theta10.bGrid, theta10.q(fixYs(3),:), theta10Style, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ High');
plot(thetad100.bGrid, thetad100.q(fixYs(1),:), thetad100Style, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ Low');
plot(thetad100.bGrid, thetad100.q(fixYs(2),:), thetad100Style, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ Med');
plot(thetad100.bGrid, thetad100.q(fixYs(3),:), thetad100Style, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ High');
xlabel('$B''$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$q(B'',y)$', 'Interpreter', 'latex', 'FontSize', 12);
% IMPORTANT: Remove legend call from here.
% legend('Location', 'best', 'Interpreter', 'latex');
grid on;
xlim([0 0.7]); % 规定横轴为从0 到 0.7
formatFigure();

%% 图4：利差对比 (Spread)
figure(4);
set(gcf, 'Position', [250, 250, 800, 600]);
sp_baseline = (baseline.delta + baseline.rf) * (1 ./ baseline.q(fixYs, :) - 1);
sp_theta10 = (theta10.delta + theta10.rf) * (1 ./ theta10.q(fixYs, :) - 1);
sp_thetad100 = (thetad100.delta + thetad100.rf) * (1 ./ thetad100.q(fixYs, :) - 1);
hold on;
plot(baseline.bGrid, (1+sp_baseline(1,:)).^4-1, baselineStyle, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline Low');
plot(baseline.bGrid, (1+sp_baseline(2,:)).^4-1, baselineStyle, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline Med');
plot(baseline.bGrid, (1+sp_baseline(3,:)).^4-1, baselineStyle, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline High');
plot(theta10.bGrid, (1+sp_theta10(1,:)).^4-1, theta10Style, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ Low');
plot(theta10.bGrid, (1+sp_theta10(2,:)).^4-1, theta10Style, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ Med');
plot(theta10.bGrid, (1+sp_theta10(3,:)).^4-1, theta10Style, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ High');
plot(thetad100.bGrid, (1+sp_thetad100(1,:)).^4-1, thetad100Style, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ Low');
plot(thetad100.bGrid, (1+sp_thetad100(2,:)).^4-1, thetad100Style, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ Med');
plot(thetad100.bGrid, (1+sp_thetad100(3,:)).^4-1, thetad100Style, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ High');
xlabel('$B$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Spread', 'Interpreter', 'latex', 'FontSize', 12);
% IMPORTANT: Remove legend call from here.
% legend('Location', 'best', 'Interpreter', 'latex');
grid on;
xlim([0 0.7]); % 规定横轴为从0 到 0.7
ylim([0 0.5]); % 规定利差图的纵轴从 0 到 0.2
formatFigure();

%% 图5：B* 交叉点示意图 (Pivoting Bond Price Schedules)
% 为了展示"doom softening"效应，使用一个正常的产出水平
figure(5);
set(gcf, 'Position', [300, 300, 800, 600]);
hold on;

% 使用中等产出水平（正常状态）
normalY = fixYs(2);

% 绘制三种theta值下的债券价格曲线
plot(baseline.bGrid, baseline.q(normalY,:), baselineStyle, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', '$\theta = 1$ (Baseline)');
plot(theta10.bGrid, theta10.q(normalY,:), theta10Style, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', '$\theta = 10$ (Medium)');
plot(thetad100.bGrid, thetad100.q(normalY,:), thetad100Style, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', '$\theta = 100$ (High)');

% 找到交叉点 B*
% 寻找 baseline 和 thetad100 曲线的交叉点
q_diff = baseline.q(normalY,:) - thetad100.q(normalY,:);
% 找到符号变化的位置
sign_changes = find(diff(sign(q_diff)) ~= 0);
if ~isempty(sign_changes)
    % 使用第一个交叉点
    cross_idx = sign_changes(1);
    if cross_idx < length(baseline.bGrid)
        % 线性插值找到更精确的交叉点
        x1 = baseline.bGrid(cross_idx);
        x2 = baseline.bGrid(cross_idx+1);
        y1 = q_diff(cross_idx);
        y2 = q_diff(cross_idx+1);
        B_star = x1 - y1 * (x2 - x1) / (y2 - y1);
        q_star = interp1(baseline.bGrid, baseline.q(normalY,:), B_star);

        % 标记交叉点
        plot(B_star, q_star, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red', 'DisplayName', '$B^*(y)$');

        % 添加垂直线标记交叉点
        plot([B_star B_star], [0 q_star], 'r--', 'LineWidth', 1, 'HandleVisibility', 'off');

        % 在图上标注区域
        text(B_star/2, 0.85, 'PRO Premium', 'FontSize', 14, 'Interpreter', 'latex', ...
            'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'EdgeColor', 'black');
        text((B_star + 0.7)/2, 0.85, 'Doom Softening', 'FontSize', 14, 'Interpreter', 'latex', ...
            'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'EdgeColor', 'black');
    end
end

xlabel('$B''$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$q(B'',y)$', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
xlim([0 0.7]);
ylim([0 1]);
formatFigure();

%% 保存图像
fprintf('\nSaving figures 1-5...\n');
for i = 1:5
    fig = figure(i);
    set(findall(fig, '-property', 'FontName'), 'FontName', 'Times New Roman');
    % formatFigure() is already called within each figure block
    set(fig, 'PaperPositionMode', 'auto');
    fig_pos = get(fig, 'PaperPosition');
    set(fig, 'PaperSize', [fig_pos(3) fig_pos(4)]);
    filename = sprintf('comparison_figure_%d.pdf', i);
    print(fig, filename, '-dpdf', '-r300');
end

%% 图10：债务阈值 B*(y) 的单调性 (Monotonicity of Debt Threshold)
% 展示Proposition 2: B*(y) 随 y 单调递增
figure(10);
set(gcf, 'Position', [350, 350, 800, 600]);
hold on;

% 计算不同收入水平下的交叉点 B*(y)
n_y_points = 50; % 选择20个收入水平点
y_indices = round(linspace(10, baseline.ySz-10, n_y_points)); % 避免极端值
B_star_baseline_theta10 = zeros(n_y_points, 1);
B_star_baseline_thetad100 = zeros(n_y_points, 1);
y_values = zeros(n_y_points, 1);

for i = 1:n_y_points
    yi = y_indices(i);
    y_values(i) = baseline.yGrid(yi);

    % 计算 baseline 和 theta10 的交叉点
    q_diff_10 = baseline.q(yi,:) - theta10.q(yi,:);
    sign_changes_10 = find(diff(sign(q_diff_10)) ~= 0);
    if ~isempty(sign_changes_10)
        cross_idx = sign_changes_10(1);
        if cross_idx < length(baseline.bGrid)
            % 线性插值找到交叉点
            x1 = baseline.bGrid(cross_idx);
            x2 = baseline.bGrid(cross_idx+1);
            y1 = q_diff_10(cross_idx);
            y2 = q_diff_10(cross_idx+1);
            B_star_baseline_theta10(i) = x1 - y1 * (x2 - x1) / (y2 - y1);
        else
            B_star_baseline_theta10(i) = NaN;
        end
    else
        B_star_baseline_theta10(i) = NaN;
    end

    % 计算 baseline 和 thetad100 的交叉点
    q_diff_100 = baseline.q(yi,:) - thetad100.q(yi,:);
    sign_changes_100 = find(diff(sign(q_diff_100)) ~= 0);
    if ~isempty(sign_changes_100)
        cross_idx = sign_changes_100(1);
        if cross_idx < length(baseline.bGrid)
            % 线性插值找到交叉点
            x1 = baseline.bGrid(cross_idx);
            x2 = baseline.bGrid(cross_idx+1);
            y1 = q_diff_100(cross_idx);
            y2 = q_diff_100(cross_idx+1);
            B_star_baseline_thetad100(i) = x1 - y1 * (x2 - x1) / (y2 - y1);
        else
            B_star_baseline_thetad100(i) = NaN;
        end
    else
        B_star_baseline_thetad100(i) = NaN;
    end
end

% 移除NaN值
valid_10 = ~isnan(B_star_baseline_theta10);
valid_100 = ~isnan(B_star_baseline_thetad100);

% 绘制 B*(y) 曲线
plot(y_values(valid_10), B_star_baseline_theta10(valid_10), theta10Style, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', '$B^*_{1,10}(y)$');
plot(y_values(valid_100), B_star_baseline_thetad100(valid_100), thetad100Style, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', '$B^*_{1,100}(y)$');

% 添加标记点以强调单调性
% scatter(y_values(valid_10), B_star_baseline_theta10(valid_10), 50, colors(2,:), 'filled', 'HandleVisibility', 'off');
% scatter(y_values(valid_100), B_star_baseline_thetad100(valid_100), 50, colors(3,:), 'filled', 'HandleVisibility', 'off');

% 添加趋势线以突出单调性
if sum(valid_10) > 1
    p_10 = polyfit(y_values(valid_10), B_star_baseline_theta10(valid_10), 1);
    y_trend_10 = linspace(min(y_values(valid_10)), max(y_values(valid_10)), 100);
    B_trend_10 = polyval(p_10, y_trend_10);
    plot(y_trend_10, B_trend_10, ':', 'Color', colors(2,:), 'LineWidth', 1, 'HandleVisibility', 'off');
end

if sum(valid_100) > 1
    p_100 = polyfit(y_values(valid_100), B_star_baseline_thetad100(valid_100), 1);
    y_trend_100 = linspace(min(y_values(valid_100)), max(y_values(valid_100)), 100);
    B_trend_100 = polyval(p_100, y_trend_100);
    plot(y_trend_100, B_trend_100, ':', 'Color', colors(3,:), 'LineWidth', 1, 'HandleVisibility', 'off');
end

xlabel('$y$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$B^*(y)$', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
formatFigure();

hold off;

%% 新增对比图
fprintf('\nCreating additional comparison plots...\n');

% 图6：利差分布直方图
figure(6);
set(gcf, 'Position', [100, 100, 800, 600]);
hold on;
histogram(baseline.spSim(baseline.valid)*4, 50, 'Normalization', 'probability', 'DisplayName', 'Baseline', 'FaceAlpha', 0.5, 'FaceColor', colors(1,:), 'EdgeColor', 'none');
histogram(theta10.spSim(theta10.valid)*4, 50, 'Normalization', 'probability', 'DisplayName', 'Med $\theta$', 'FaceAlpha', 0.5, 'FaceColor', colors(2,:), 'EdgeColor', 'none');
histogram(thetad100.spSim(thetad100.valid)*4, 50, 'Normalization', 'probability', 'DisplayName', 'High $\theta$', 'FaceAlpha', 0.5, 'FaceColor', colors(3,:), 'EdgeColor', 'none');
hold off;
xlabel('Spread', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Frequency', 'Interpreter', 'latex', 'FontSize', 12);
% IMPORTANT: Remove legend call from here.
% legend('Location', 'best', 'Interpreter', 'latex');
grid on;
formatFigure();

% 图7：债务分布直方图
figure(7);
set(gcf, 'Position', [150, 150, 800, 600]);
% 计算年化的债务/GDP比率序列
debt_gdp_baseline = baseline.bSim(baseline.valid) ./ (baseline.gdpSim(baseline.valid) * 4);
debt_gdp_theta10  = theta10.bSim(theta10.valid) ./ (theta10.gdpSim(theta10.valid) * 4);
debt_gdp_thetad100 = thetad100.bSim(thetad100.valid) ./ (thetad100.gdpSim(thetad100.valid) * 4);

hold on;
% 绘制比率序列的直方图
histogram(debt_gdp_baseline, 50, 'Normalization', 'probability', 'DisplayName', 'Baseline', 'FaceAlpha', 0.5, 'FaceColor', colors(1,:), 'EdgeColor', 'none');
histogram(debt_gdp_theta10, 50, 'Normalization', 'probability', 'DisplayName', 'Med $\theta$', 'FaceAlpha', 0.5, 'FaceColor', colors(2,:), 'EdgeColor', 'none');
histogram(debt_gdp_thetad100, 50, 'Normalization', 'probability', 'DisplayName', 'High $\theta$', 'FaceAlpha', 0.5, 'FaceColor', colors(3,:), 'EdgeColor', 'none');
hold off;
xlabel('Debt/GDP', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Frequency', 'Interpreter', 'latex', 'FontSize', 12);
% IMPORTANT: Remove legend call from here.
% legend('Location', 'best', 'Interpreter', 'latex');
grid on;
formatFigure();

% 图8：债务选择策略函数（Choice Probability）
figure(8);
set(gcf, 'Position', [200, 200, 800, 600]);
hold on;
b_idx = round(baseline.bSz * 0.6); % 使用动态索引，大约在60%位置

% Baseline
plot(baseline.bGrid, squeeze(baseline.bPol(fixYs(1), b_idx, :)), baselineStyle, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline Low');
plot(baseline.bGrid, squeeze(baseline.bPol(fixYs(2), b_idx, :)), baselineStyle, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline Med');
plot(baseline.bGrid, squeeze(baseline.bPol(fixYs(3), b_idx, :)), baselineStyle, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline High');
% Med theta
plot(theta10.bGrid, squeeze(theta10.bPol(fixYs(1), b_idx, :)), theta10Style, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ Low');
plot(theta10.bGrid, squeeze(theta10.bPol(fixYs(2), b_idx, :)), theta10Style, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ Med');
plot(theta10.bGrid, squeeze(theta10.bPol(fixYs(3), b_idx, :)), theta10Style, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ High');
% High theta
plot(thetad100.bGrid, squeeze(thetad100.bPol(fixYs(1), b_idx, :)), thetad100Style, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ Low');
plot(thetad100.bGrid, squeeze(thetad100.bPol(fixYs(2), b_idx, :)), thetad100Style, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ Med');
plot(thetad100.bGrid, squeeze(thetad100.bPol(fixYs(3), b_idx, :)), thetad100Style, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ High');
hold off;
xlabel('$B''$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Probability', 'Interpreter', 'latex', 'FontSize', 12);
% IMPORTANT: Get handle for title only if necessary, otherwise handle legend in formatFigure()
% lgd = legend('Location', 'best', 'Interpreter', 'latex');
% title(lgd, sprintf('For Current B = %.2f', baseline.bGrid(b_idx)));
grid on;
xlim([0.18 0.3]);
formatFigure(true, sprintf('For Current B = %.2f', baseline.bGrid(b_idx))); % Pass title to formatFigure

% 图9：预期下期债务 (Expected Next Period Debt)
figure(9);
set(gcf, 'Position', [250, 250, 800, 600]);
hold on;
% Baseline
plot(baseline.bGrid, baseline.EbPr(fixYs(1),:), baselineStyle, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline Low');
plot(baseline.bGrid, baseline.EbPr(fixYs(2),:), baselineStyle, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline Med');
plot(baseline.bGrid, baseline.EbPr(fixYs(3),:), baselineStyle, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline High');
% Med theta
plot(theta10.bGrid, theta10.EbPr(fixYs(1),:), theta10Style, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ Low');
plot(theta10.bGrid, theta10.EbPr(fixYs(2),:), theta10Style, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ Med');
plot(theta10.bGrid, theta10.EbPr(fixYs(3),:), theta10Style, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ High');
% High theta
plot(thetad100.bGrid, thetad100.EbPr(fixYs(1),:), thetad100Style, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ Low');
plot(thetad100.bGrid, thetad100.EbPr(fixYs(2),:), thetad100Style, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ Med');
plot(thetad100.bGrid, thetad100.EbPr(fixYs(3),:), thetad100Style, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ High');
plot(baseline.bGrid, baseline.bGrid, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
hold off;
xlabel('$B$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\mathrm{E}[B''|y,B]$', 'Interpreter', 'latex', 'FontSize', 12);
% IMPORTANT: Remove legend call from here.
% legend('Location', 'best', 'Interpreter', 'latex');
grid on;
xlim([0 0.35])
formatFigure();

% 保存新增的图像 (6-10)
fprintf('\nSaving figures 6-10...\n');
for i = 6:10
    fprintf('  Saving figure %d...\n', i);
    fig = figure(i);
    set(findall(fig, '-property', 'FontName'), 'FontName', 'Times New Roman');
    % formatFigure() is already called within each figure block
    set(fig, 'PaperPositionMode', 'auto');
    fig_pos = get(fig, 'PaperPosition');
    set(fig, 'PaperSize', [fig_pos(3) fig_pos(4)]);
    filename = sprintf('comparison_figure_%d.pdf', i);
    print(fig, filename, '-dpdf', '-r300');
    fprintf('  Figure %d saved successfully.\n', i);
end

fprintf('\nAll additional comparison plots saved.\n');


%% 图11 & 12：去杠杆路径对比 (Deleveraging Dynamics)
fprintf('\nCreating deleveraging dynamics plots...\n');

% --- 1. 模拟设置 ---
T_sim = 50;         % 模拟的总期数（季度）
B0 = 0.35;           % 初始债务水平
fprintf('Using y_idx = %d, corresponding to y = %.4f\n', fixYs(2), baseline.yGrid(fixYs(2)));

% --- 2. 模拟三个经济体的去杠杆过程 ---
% 调用辅助函数进行模拟
deleveraging_baseline = simulate_path(baseline, B0, fixYs(2), T_sim);
deleveraging_theta10  = simulate_path(theta10, B0, fixYs(2), T_sim);
deleveraging_thetad100 = simulate_path(thetad100, B0, fixYs(2), T_sim);

deleveraging_baseline_low  = simulate_path(baseline, B0, fixYs(1), T_sim);
deleveraging_theta10_low  = simulate_path(theta10, B0, fixYs(1), T_sim);
deleveraging_thetad100_low = simulate_path(thetad100, B0, fixYs(1), T_sim);

deleveraging_baseline_high  = simulate_path(baseline, B0, fixYs(3), T_sim);
deleveraging_theta10_high  = simulate_path(theta10, B0, fixYs(3), T_sim);
deleveraging_thetad100_high = simulate_path(thetad100, B0, fixYs(3), T_sim);

% 时间轴（以年为单位）
time_years = (0:T_sim-1) / 4;

%% 图11：债务路径对比 (Debt Path)
figure(11);
set(gcf, 'Position', [100, 100, 800, 600]);
hold on;
plot(time_years, deleveraging_baseline.B_path, baselineStyle, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline Med');
plot(time_years, deleveraging_theta10.B_path, theta10Style, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ Med');
plot(time_years, deleveraging_thetad100.B_path, thetad100Style, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ Med');
plot(time_years, deleveraging_baseline_low.B_path, baselineStyle, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline Low');
plot(time_years, deleveraging_theta10_low.B_path, theta10Style, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ Low');
plot(time_years, deleveraging_thetad100_low.B_path, thetad100Style, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ Low');
plot(time_years, deleveraging_baseline_high.B_path, baselineStyle, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline High');
plot(time_years, deleveraging_theta10_high.B_path, theta10Style, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ High');
plot(time_years, deleveraging_thetad100_high.B_path, thetad100Style, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ High');
hold off;

xlabel('Quarters', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$B_t$', 'Interpreter', 'latex', 'FontSize', 12);

grid on;
xlim([0 6]);
legend('Location', 'best', 'Interpreter', 'latex');
formatFigure(); % 使用一个简化的格式化函数

%% 图12：消费路径对比 (Consumption Path)
figure(12);
set(gcf, 'Position', [150, 150, 800, 600]);
hold on;
plot(time_years, deleveraging_baseline.C_path, baselineStyle, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline Med');
plot(time_years, deleveraging_theta10.C_path, theta10Style, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ Med');
plot(time_years, deleveraging_thetad100.C_path, thetad100Style, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ Med');
plot(time_years, deleveraging_baseline_low.C_path, baselineStyle, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline Low');
plot(time_years, deleveraging_theta10_low.C_path, theta10Style, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ Low');
plot(time_years, deleveraging_thetad100_low.C_path, thetad100Style, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ Low');
plot(time_years, deleveraging_baseline_high.C_path, baselineStyle, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline High');
plot(time_years, deleveraging_theta10_high.C_path, theta10Style, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ High');
plot(time_years, deleveraging_thetad100_high.C_path, thetad100Style, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ High');
hold off;

xlabel('Quarters', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$c_t$', 'Interpreter', 'latex', 'FontSize', 12);

grid on;
xlim([0 6]);
legend('Location', 'best', 'Interpreter', 'latex');
formatFigure(); % 使用一个简化的格式化函数

figure(21);
set(gcf, 'Position', [200, 200, 800, 600]);
hold on;
% 绘制中等产出水平下的利差路径
plot(time_years, deleveraging_baseline.Sp_path, baselineStyle, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline Med');
plot(time_years, deleveraging_theta10.Sp_path, theta10Style, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ Med');
plot(time_years, deleveraging_thetad100.Sp_path, thetad100Style, 'Color', colors(2,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ Med');
% (可选) 绘制低产出水平下的利差路径
plot(time_years, deleveraging_baseline_low.Sp_path, baselineStyle, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline Low');
plot(time_years, deleveraging_theta10_low.Sp_path, theta10Style, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ Low');
plot(time_years, deleveraging_thetad100_low.Sp_path, thetad100Style, 'Color', colors(1,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ Low');
% (可选) 绘制高产出水平下的利差路径
plot(time_years, deleveraging_baseline_high.Sp_path, baselineStyle, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'Baseline High');
plot(time_years, deleveraging_theta10_high.Sp_path, theta10Style, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'Med $\theta$ High');
plot(time_years, deleveraging_thetad100_high.Sp_path, thetad100Style, 'Color', colors(3,:), 'LineWidth', lineWidth, 'DisplayName', 'High $\theta$ High');
hold off;

xlabel('Quarters', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Annualized Spread', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
xlim([0 6]);
legend('Location', 'best', 'Interpreter', 'latex');
formatFigure(); % 使用您已有的格式化函数

%% 保存图像 (11-12)
fprintf('\nSaving deleveraging dynamics figures (11,12,21)...\n');
for i = [11,12,21]
    fig = figure(i);
    set(findall(fig, '-property', 'FontName'), 'FontName', 'Times New Roman');
    set(fig, 'PaperPositionMode', 'auto');
    fig_pos = get(fig, 'PaperPosition');
    set(fig, 'PaperSize', [fig_pos(3) fig_pos(4)]);
    filename = sprintf('comparison_figure_%d.pdf', i);
    print(fig, filename, '-dpdf', '-r300');
end
fprintf('Figures 11, 12, 21 saved successfully.\n');


%% 图13-16：对产出冲击的脉冲响应 (Impulse Response Functions)
fprintf('\nCreating Impulse Response Function plots...\n');

% --- 1. 模拟设置 ---
T_irf = 24;          % IRF的模拟期数 (季度), 10年
y_mean = 1.0;        % 产出均值
shock_size = 0.03;   % 产出冲击的大小 (3%)
y_shock_level = y_mean * (1 + shock_size); % 冲击后的产出水平

% 找到最接近均值和冲击后产出的网格点索引
[~, y_idx_mean] = min(abs(baseline.yGrid - y_mean));
[~, y_idx_shock] = min(abs(baseline.yGrid - y_shock_level));

fprintf('Using y_idx_mean = %d (y=%.4f) and y_idx_shock = %d (y=%.4f)\n', ...
    y_idx_mean, baseline.yGrid(y_idx_mean), y_idx_shock, baseline.yGrid(y_idx_shock));

% 构建产出冲击路径 (t=1时有冲击，之后恢复正常)
y_path_idx = ones(T_irf, 1) * y_idx_mean;
y_path_idx(1) = y_idx_shock;

% --- 2. 为每个经济体寻找不动点并模拟IRF ---
% 调用辅助函数
irf_baseline = find_ss_and_simulate_irf(baseline, y_path_idx);
irf_theta10  = find_ss_and_simulate_irf(theta10, y_path_idx);
irf_thetad100 = find_ss_and_simulate_irf(thetad100, y_path_idx);

% 时间轴（季度）
time_quarters = 0:T_irf-1;

% 定义一个函数来绘制每张图，避免代码重复
plot_irf(13, time_quarters, '$y_t$', ...
    {irf_baseline.y_path, irf_theta10.y_path, irf_thetad100.y_path}, ...
    'Deviation from ss. $(y_t, \%)$', baselineStyle, theta10Style, thetad100Style, colors, lineWidth);

plot_irf(14, time_quarters, '$B_t$', ...
    {irf_baseline.B_path, irf_theta10.B_path, irf_thetad100.B_path}, ...
    'Deviation from ss. $(B_t, \%)$', baselineStyle, theta10Style, thetad100Style, colors, lineWidth);

plot_irf(15, time_quarters, '$c_t$', ...
    {irf_baseline.C_path, irf_theta10.C_path, irf_thetad100.C_path}, ...
    'Deviation from ss. $(c_t, \%)$', baselineStyle, theta10Style, thetad100Style, colors, lineWidth);

% 注意：利差是年化的
plot_irf(16, time_quarters, 'Spread', ...
    {irf_baseline.Sp_path, irf_theta10.Sp_path, irf_thetad100.Sp_path}, ...
    'Deviation from ss. (Spread, bps)', baselineStyle, theta10Style, thetad100Style, colors, lineWidth);

%% 保存图像 (13-16)
fprintf('\nSaving Transitory IRF figures (13-16)...\n');
for i = 13:16
    fig = figure(i);
    set(findall(fig, '-property', 'FontName'), 'FontName', 'Times New Roman');
    set(fig, 'PaperPositionMode', 'auto');
    fig_pos = get(fig, 'PaperPosition');
    set(fig, 'PaperSize', [fig_pos(3) fig_pos(4)]);
    filename = sprintf('comparison_figure_%d.pdf', i);
    print(fig, filename, '-dpdf', '-r300');
end
fprintf('Figures 13-16 saved successfully.\n');

%% 17-20：对持续性产出冲击的脉冲响应 (IRF to Persistent Shock)
fprintf('\nCreating IRF plots for a persistent shock...\n');

% --- 1. 模拟设置 ---
T_irf_persistent = 24; % 模拟期数，对于持续性冲击需要更长的时间
rho_y = 0.8;          % AR(1) 冲击的自回归系数
initial_shock_size = 0.03; % 初始冲击大小 (3%)
y_mean = 1.0;          % 产出均值

% --- 2. 生成AR(1)冲击路径 ---
% 初始化冲击序列 (对数偏离)
log_y_shock_path = zeros(T_irf_persistent, 1);
log_y_shock_path(1) = log(1 + initial_shock_size); % 初始对数冲击

% 迭代生成AR(1)路径
for t = 2:T_irf_persistent
    log_y_shock_path(t) = rho_y * log_y_shock_path(t-1);
end

% 将对数冲击转换为产出水平路径
y_path_persistent = y_mean * exp(log_y_shock_path);

% 将连续的产出路径映射到离散的网格索引上
y_path_idx_persistent = zeros(T_irf_persistent, 1);
for t = 1:T_irf_persistent
    [~, y_path_idx_persistent(t)] = min(abs(baseline.yGrid - y_path_persistent(t)));
end

fprintf('Generated persistent shock path with rho = %.2f\n', rho_y);

% --- 3. 为每个经济体寻找不动点并模拟IRF ---
% 复用之前的辅助函数，只需传入新的冲击路径
irf_p_baseline = find_ss_and_simulate_irf(baseline, y_path_idx_persistent);
irf_p_theta10  = find_ss_and_simulate_irf(theta10, y_path_idx_persistent);
irf_p_thetad100 = find_ss_and_simulate_irf(thetad100, y_path_idx_persistent);

% --- 4. 绘图 ---
% 定义线条样式和颜色


% 时间轴（季度）
time_quarters_p = 0:T_irf_persistent-1;

% 定义一个函数来绘制每张图，避免代码重复
% (复用之前的 plot_irf 函数)
plot_irf(17, time_quarters_p, 'Output ($y_t$)', ...
    {irf_p_baseline.y_path, irf_p_theta10.y_path, irf_p_thetad100.y_path}, ...
    'Deviation from ss. $(y_t, \%)$', baselineStyle, theta10Style, thetad100Style, colors, lineWidth);

plot_irf(18, time_quarters_p, 'Debt Stock ($B_t$)', ...
    {irf_p_baseline.B_path, irf_p_theta10.B_path, irf_p_thetad100.B_path}, ...
    'Deviation from ss. $(B_t, \%)$', baselineStyle, theta10Style, thetad100Style, colors, lineWidth);

plot_irf(19, time_quarters_p, 'Consumption ($c_t$)', ...
    {irf_p_baseline.C_path, irf_p_theta10.C_path, irf_p_thetad100.C_path}, ...
    'Deviation from ss. $(c_t, \%)$', baselineStyle, theta10Style, thetad100Style, colors, lineWidth);

% 注意：利差是年化的
plot_irf(20, time_quarters_p, 'Spread', ...
    {irf_p_baseline.Sp_path, irf_p_theta10.Sp_path, irf_p_thetad100.Sp_path}, ...
    'Deviation from ss. (Spread, bps)', baselineStyle, theta10Style, thetad100Style, colors, lineWidth);

%% 保存图像 (17-20)
fprintf('\nSaving Persistent IRF figures (17-20)...\n');
for i = 17:20
    fig = figure(i);
    set(findall(fig, '-property', 'FontName'), 'FontName', 'Times New Roman');
    set(fig, 'PaperPositionMode', 'auto');
    fig_pos = get(fig, 'PaperPosition');
    set(fig, 'PaperSize', [fig_pos(3) fig_pos(4)]);
    filename = sprintf('comparison_figure_%d.pdf', i);
    print(fig, filename, '-dpdf', '-r300');
end
fprintf('Figures 17-20 saved successfully.\n');








%% --- 辅助函数 ---

function irf_data = find_ss_and_simulate_irf(model, y_path_idx)
    % 1. 寻找不动点 B_ss
    y_idx_mean = y_path_idx(end); % 稳态时的y索引

    % 定义函数 g(B) = E[B'|y_mean, B] - B
    expected_b_prime_func = @(b_idx) sum(squeeze(model.bPol(y_idx_mean, b_idx, :)) .* model.bGrid);
    g_func = @(b_idx) expected_b_prime_func(b_idx) - model.bGrid(b_idx);

    % 在整个债务网格上计算 g(B)
    g_values = arrayfun(g_func, 1:model.bSz);

    % 找到 g(B) 符号变化的位置，即不动点所在
    sign_changes = find(diff(sign(g_values)) ~= 0);
    if isempty(sign_changes)
        % 如果没有交叉点，可能在边界，取g值最接近0的点
        [~, ss_idx] = min(abs(g_values));
        warning('No steady-state crossing found. Using point closest to zero.');
    else
        % 线性插值找到更精确的不动点
        idx1 = sign_changes(1);
        idx2 = idx1 + 1;
        B_ss = interp1(g_values([idx1, idx2]), model.bGrid([idx1, idx2]), 0);
        [~, ss_idx] = min(abs(model.bGrid - B_ss)); % 找到最近的网格点
    end

    B_ss = model.bGrid(ss_idx);
    fprintf('Found steady state debt B_ss = %.4f for model.\n', B_ss);

    % 2. 模拟IRF路径
    T = length(y_path_idx);
    irf_data = struct('y_path', zeros(T,1), 'B_path', zeros(T,1), ...
        'C_path', zeros(T,1), 'Sp_path', zeros(T,1));

    % 设置初始状态
    b_idx_current = ss_idx;

    for t = 1:T
        % 获取当前状态
        B_current = model.bGrid(b_idx_current);
        y_idx_current = y_path_idx(t);
        y_current = model.yGrid(y_idx_current);

        % 确定下一期债务 B' (取期望)
        bPol_current = squeeze(model.bPol(y_idx_current, b_idx_current, :));
        B_next = sum(bPol_current .* model.bGrid);
        [~, b_idx_next] = min(abs(model.bGrid - B_next));

        % 获取价格 q(B', y)
        q_next = model.q(y_idx_current, b_idx_next);

        % 计算当期消费 C_t 和 利差 Sp_t
        kappa = model.delta + model.rf;
        C_t = y_current - kappa * B_current + q_next * (B_next - (1 - model.delta) * B_current);
        Sp_t = (kappa * (1/q_next - 1)) * 4; % 算术年化

        % 存储结果
        irf_data.y_path(t) = y_current;
        irf_data.B_path(t) = B_current;
        irf_data.C_path(t) = C_t;
        irf_data.Sp_path(t) = Sp_t;

        % 更新状态
        b_idx_current = b_idx_next;
    end

    % 3. 计算稳态值并报告偏离
    y_ss = model.yGrid(y_idx_mean);
    C_ss = irf_data.C_path(end); % 假设T期后已回到稳态
    Sp_ss = irf_data.Sp_path(end);

    irf_data.y_path = (irf_data.y_path - y_ss) / y_ss; % 产出：百分比偏离
    irf_data.C_path = (irf_data.C_path - C_ss) / C_ss * 100; % 消费：百分比偏离
    irf_data.Sp_path = (irf_data.Sp_path - Sp_ss) * 10000; % 利差：基点(bps)偏离
    % 债务报告水平值，因为它本身就是我们关心的状态变量
    irf_data.B_path = (irf_data.B_path - B_ss) / B_ss * 100;
end

function plot_irf(fig_num, time, y_label, data_cell, y_axis_label, ls1, ls2, ls3, colors, lw)
    figure(fig_num);
    set(gcf, 'Position', [100 + 50*(fig_num-13), 100 + 50*(fig_num-13), 800, 600]);
    hold on;
    plot(time, data_cell{1}, ls1, 'Color', colors(2,:), 'LineWidth', lw, 'DisplayName', 'Baseline');
    plot(time, data_cell{2}, ls2, 'Color', colors(2,:), 'LineWidth', lw, 'DisplayName', 'Med $\theta$');
    plot(time, data_cell{3}, ls3, 'Color', colors(2,:), 'LineWidth', lw, 'DisplayName', 'High $\theta$');
    plot(time, zeros(size(time)), 'k-', 'LineWidth', 1, 'HandleVisibility', 'off'); % 零线
    hold off;

    xlabel('Quarters', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel(y_axis_label, 'Interpreter', 'latex', 'FontSize', 12);
    grid on;
    xlim([0, time(end)]);
    legend('Location', 'best', 'Interpreter', 'latex');
    formatFigure();
end




%% --- 辅助函数 ---

function path = simulate_path(model, B0, y_idx, T)
    % 该函数模拟给定初始债务和固定产出下的路径，并正确记录所有变量

    % 初始化结果存储
    path.B_path = zeros(T, 1);
    path.C_path = zeros(T, 1);
    path.Sp_path = zeros(T, 1); % 新增利差路径的初始化

    % 获取模型参数和网格
    bGrid = model.bGrid;
    y = model.yGrid(y_idx);
    kappa = model.delta + model.rf;

    % 找到最接近初始债务 B0 的网格点索引
    [~, b_idx_current] = min(abs(bGrid - B0));
    path.B_path(1) = bGrid(b_idx_current);

    % --- 循环结构修正 ---
    % 循环 T 次，每次计算 t 时刻的 C 和 Sp，以及 t+1 时刻的 B
    for t = 1:T
        % 1. 获取当前状态 (y_idx, b_idx_current)
        B_current = path.B_path(t);
        [~, b_idx_current] = min(abs(bGrid - B_current));

        % 2. 确定下一期债务 B' (取期望)
        bPol_current = squeeze(model.bPol(y_idx, b_idx_current, :));
        B_next = sum(bPol_current .* bGrid);
        [~, b_idx_next] = min(abs(bGrid - B_next));

        % 3. 获取价格 q(B', y) 并计算当期利差
        q_next = model.q(y_idx, b_idx_next);
        Sp_t_quarterly = (kappa * (1/q_next - 1));
        % 与图4保持一致，使用几何年化
        path.Sp_path(t) = (1 + Sp_t_quarterly)^4 - 1;

        % 4. 计算当期消费 C_t
        path.C_path(t) = y - kappa * B_current + q_next * (B_next - (1 - model.delta) * B_current);

        % 5. 存储下一期债务（如果不是最后一次循环）
        if t < T
            path.B_path(t+1) = B_next;
        end
    end
end

% --- MODIFIED formatFigure function ---
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