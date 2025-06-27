using Interpolations, Plots
using Random
Random.seed!(123)

function interp(x, y)
    itp1d = interpolate((x,), y, Gridded(Linear()))
    eitp1d = extrapolate(itp1d, Line())
    return eitp1d
end

# 一维数据（线性插值+线性外推）
x = 1.0:5.0
y = 1 * x .+ rand(5)
eitp1d = interp(x, y)

eitp_array = [eitp1d(i) for i = 1.0:10.0]

plot(
    x,
    y,
    label = "data",
    xlabel = "x",
    ylabel = "y",
    title = "interp",
    xlims = (0, 10),
    ylims = (0, 10),
)
plot!(1.0:10.0, eitp_array, label = "interp")
savefig("Scripts/interp.pdf")

# eitp1d(11.0)
# eitp1d(100.0)

# 二维数据（线性插值+线性外推）
xs = 1.0:10.0
ys = 1.0:10.0
A = 1 * xs .+ 1 * ys' .+ rand(10, 10)
function interp2d(xs, ys, A)
    itp2d = interpolate((xs, ys), A, Gridded(Linear()))
    eitp2d = extrapolate(itp2d, Line())
    return eitp2d
end

eitp2d = interp2d(xs, ys, A)
println("二维插值点(5.5, 7.2): ", eitp2d(5.5, 7.2))
println("二维外推点(0.0, 0.0): ", eitp2d(0.0, 0.0))
println("二维外推点(11.0, 11.0): ", eitp2d(11.0, 11.0))

# 创建网格点用于3D绘图
x_grid = range(0, 10, length = 100)
y_grid = range(0, 10, length = 100)
z_grid = [eitp2d(x, y) for x in x_grid, y in y_grid]

# 绘制3D表面图
surface(
    x_grid,
    y_grid,
    z_grid,
    xlabel = "x",
    ylabel = "y",
    zlabel = "z",
    title = "2D Interpolation Surface",
    camera = (30, 30),
)  # 设置视角
savefig("Scripts/interp2d_3d.pdf")
