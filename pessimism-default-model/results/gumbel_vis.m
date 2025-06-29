% --- 生成Gumbel分布图的MATLAB代码 ---
% 定义三种颜色，并将其从RGB 0-255范围转换到0-1范围
colors = [188 228 183;     % 用于PDF的颜色
    064 171 092;     % (备用颜色)
    000 078 045];    % 用于CDF的颜色
colors = colors / 255;

% 设置工作目录 (可选)
cd(pwd);

% --- Gumbel分布的参数 ---
eta = 5e-4; % 尺度参数 (Scale parameter)
euler_c = 0.5772156649; % 欧拉-马斯切罗尼常数 (Euler-Mascheroni constant)
mu_L = -eta * euler_c; % 位置参数, 用于确保均值为零 (Location parameter)

% 定义 taste shock 变量的范围
% 范围与 eta 成比例
epsilon = linspace(-4*eta, 8*eta, 400);

% 计算概率密度函数 (PDF)
pdf_gumbel = (1/eta) * exp(-(epsilon - mu_L)/eta) .* exp(-exp(-(epsilon - mu_L)/eta));

% 计算累积分布函数 (CDF)
cdf_gumbel = exp(-exp(-(epsilon - mu_L)/eta));

% --- 绘图 ---
figure('Name', 'Gumbel Distribution');

% 激活左侧 Y 轴
yyaxis left;
% 绘制 PDF 曲线
plot(epsilon, pdf_gumbel, '-', 'LineWidth', 2, 'Color', colors(2,:));
% ylabel('概率密度函数 (PDF)'); % y轴标签可以取消注释

% *** 新增代码: 设置左侧 Y 轴颜色 ***
ax = gca; % 获取当前坐标区
ax.YColor = colors(2,:); % 将Y轴颜色设置为与PDF曲线相同的颜色

% 激活右侧 Y 轴
yyaxis right;
% 绘制 CDF 曲线
plot(epsilon, cdf_gumbel, 'LineWidth', 2, 'Color', colors(3,:));
% ylabel('累积分布函数 (CDF)'); % y轴标签可以取消注释

% *** 新增代码: 设置右侧 Y 轴颜色 ***
ax = gca; % 获取当前坐标区 (此时是右侧)
ax.YColor = colors(3,:); % 将Y轴颜色设置为与CDF曲线相同的颜色


% --- 添加图形标题和标签 ---
% title('Gumbel分布 (Gumbel Distribution)'); % 标题可以取消注释
xlabel('Taste Shock Value ($\epsilon$)', 'FontSize', 15, 'Interpreter', 'latex');
legend('PDF', 'CDF', 'Location', 'best', 'FontSize', 15, 'Interpreter', 'latex');
grid on;
box on;

% 将图形保存为 PNG 文件
saveas(gcf, 'gumbel_distribution.png');

disp('图形已保存为 gumbel_distribution.png');