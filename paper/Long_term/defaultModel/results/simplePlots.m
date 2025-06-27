% 清理工作区和命令行，并加载数据文件
clear; clc;
load results.mat;      % 加载模拟结果数据
load simulation.mat;   % 加载模拟参数

% 选择y维度中间的三个点，用于后续的绘图展示
ryIx = round(ySz/2, 0);               
fixYs = [ryIx-10,  ryIx, ryIx+10];      

%% 图1：利差(spread)分布直方图
figure;
% 绘制利差（spread）的直方图，采用概率归一化，便于观察频率分布
histogram(spSim(valid), 50, 'Normalization', 'probability');
hold on;
% 在直方图上添加一条表示平均值的垂直红虚线
xline(mean(spSim(valid)), 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
hold off;
xlim([0 0.1]);  % 设置x轴显示范围
xlabel('Spread', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);
title('Spread Distribution', 'FontSize', 14);
grid on;        % 添加网格线使图形更精致
formatFigure;   % 应用自定义的图形格式设置
% exportgraphics(gca, "spreadDist.pdf", "ContentType", "vector");

%% 图2：债务(debt)分布直方图
figure;
% 绘制债务水平的直方图，利用选择的样本数据进行概率归一化
histogram(bGrid(bSimIx(valid)), 50, 'Normalization', 'probability');
hold on;
% 添加债务平均水平的标志线，便于比较数据分布的集中趋势
xline(mean(bGrid(bSimIx(valid))), 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
hold off;
xlabel('Debt', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);
title('Debt Distribution', 'FontSize', 14);
grid on;
formatFigure;
% exportgraphics(gca, "debtDist.pdf", "ContentType", "vector");

%% 图3：q函数曲线
figure;
% 对于选定的三个y状态(fixYs)，绘制q函数随债务(bGrid)变化的曲线
plot(bGrid, squeeze(q(fixYs, :)), 'LineWidth', 2);
xlabel('Debt', 'FontSize', 12);
ylabel('q Value', 'FontSize', 12);
title('q Function Across Debt for Selected States', 'FontSize', 14);
grid on;
formatFigure;

%% 图4：利差曲线
figure;
% 计算利差: 利差(sp)由 (delta + rf) 乘以逆q值的调整计算得到
sp = (delta + rf) * (1 ./ squeeze(q(fixYs, :)) - 1);
% 通过变换 (1+sp)^4 - 1 来强化微小差异，并绘制相应的曲线
plot(bGrid, (1+sp).^4 - 1, 'LineWidth', 2);
ylim([0.0 0.1]);  % 设置y轴范围，便于观察细微变化
xlabel('Debt Next Period (B'''')', 'FontSize', 12);
ylabel('Spread', 'FontSize', 12);
title('Spread Curves for Selected States', 'FontSize', 14);
grid on;
formatFigure;
% exportgraphics(gca, "spreadCurves.pdf", "ContentType", "vector");

%% 图5：债务选择策略函数（Choice Probability）
figure;
% 绘制选定y状态下，在第150个债务水平对应的债务政策函数曲线
plot(bGrid, squeeze(bPol(fixYs, 150, :)), 'o-', 'LineWidth', 2);
hold on;
% 标出第150个债务水平的参考线
xline(bGrid(150), 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--');
hold off;
xlim([0.17, 0.24]);  % 限定x轴范围，突出显示相关区间
xlabel('Debt', 'FontSize', 12);
ylabel('Choice Probability', 'FontSize', 12);
title('Debt Policy Function', 'FontSize', 14);
grid on;
formatFigure;

%% 图6：预期下期债务 (Expected Next Period Debt) 与当前债务之间的关系
figure;
% 提取选定状态下的预期债务, 并将dPol大于0.75的部分设置为NaN以剔除异常值
tmp = squeeze(EbPr(fixYs, :));
tmp(squeeze(dPol(fixYs, :)) > 0.75) = NaN;
% 绘制预期下期债务随当前债务的变化曲线
plot(bGrid, tmp, 'LineWidth', 2);
hold on;
% 同时绘制45度参考线（实际债务等于当前债务的情况）
plot(bGrid, bGrid, '--k', 'LineWidth', 1);
hold off;
xlabel('Current Debt', 'FontSize', 12);
ylabel('Next Period Debt', 'FontSize', 12);
title('Expected Debt Transition', 'FontSize', 14);
xlim([0 0.5]);
grid on;
formatFigure;