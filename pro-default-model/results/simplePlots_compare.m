function simplePlots_compare(res_base, sim_base, res_high, sim_high)
% 以你原有 simplePlots.m 为蓝本，做如下修改：

fixYs = round(res_base.ySz/2) + [-10 0 10];

% 1. q函数曲线
figure; hold on;
for i = 1:3
    plot(res_base.bGrid, squeeze(res_base.q(fixYs(i), :)), 'Color', [1 0 0]*i/3, 'LineStyle', '--', 'LineWidth', 2.5);
    plot(res_high.bGrid, squeeze(res_high.q(fixYs(i), :)), 'Color', [1 0 0]*i/3, 'LineStyle', '-', 'LineWidth', 2.5);
end
legend({'q base 1','q base 2','q base 3','q high 1','q high 2','q high 3'});
title('q Function Comparison');
grid on;

% 2. 利差曲线
figure; hold on;
sp_base = (res_base.delta + res_base.rf) * (1 ./ squeeze(res_base.q(fixYs, :)) - 1);
sp_high = (res_high.delta + res_high.rf) * (1 ./ squeeze(res_high.q(fixYs, :)) - 1);
for i = 1:3
    plot(res_base.bGrid, (1+sp_base(i,:)).^4-1, 'LineStyle', '--', 'LineWidth', 2.5);
    plot(res_high.bGrid, (1+sp_high(i,:)).^4-1, 'LineStyle', '-', 'LineWidth', 2.5);
end
legend({'Spread base 1','Spread base 2','Spread base 3','Spread high 1','Spread high 2','Spread high 3'});
title('Spread Comparison');
grid on;

% 3. 其他图同理，按上述方式分别画六条线，虚线/实线区分，LineWidth=2.5
end
