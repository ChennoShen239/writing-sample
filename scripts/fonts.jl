using Plots

# =============================================
# GR 后端官方支持字体（内置 PostScript 名称）
# 参考：https://gr-framework.org/fonts.html
#
# Times:             "Times", "Times-Roman", "Times New Roman"
# Times Italic:      "Times-Italic"
# Times Bold:        "Times-Bold"
# Times BoldItalic:  "Times-BoldItalic"
# Helvetica:         "Helvetica", "Arial"
# Helvetica Oblique: "Helvetica-Oblique"
# Helvetica Bold:    "Helvetica-Bold"
# Helvetica BoldObl: "Helvetica-BoldOblique"
# Courier:           "Courier"
# Courier Oblique:   "Courier-Oblique"
# Courier Bold:      "Courier-Bold"
# Courier BoldObl:   "Courier-BoldOblique"
# Symbol:            "Symbol"
# Bookman:           "Bookman-Light", "Bookman-Demi"
# NewCenturySchlbk:  "NewCenturySchlbk-Roman"
# AvantGarde:        "AvantGarde-Book"
# Palatino:          "Palatino-Roman"
# ZapfChancery:      "ZapfChancery-MediumItalic"
# ZapfDingbats:      "ZapfDingbats"
#
# 你可以用 fontfamily="Times"/"Helvetica"/"Courier"/"Symbol" 等。
# =============================================

# 自动检测系统可用字体（需安装 fc-list，macOS/Linux 通用）
import Base: run, readchomp
function list_system_fonts()
    try
        fonts_raw = String(readchomp(`fc-list : family`))
        fonts = split(fonts_raw, '\n')
        fonts = [strip(split(line, ',')[1]) for line in fonts if !isempty(line)]
        fonts = unique(sort(fonts))
        return fonts
    catch e
        @warn "无法调用 fc-list，系统字体检测失败：$e"
        return String[]
    end
end

# 官方支持字体名（可自定义增减）
gr_fonts = [
    "Times",
    "Times-Roman",
    "Times New Roman",
    "Times-Italic",
    "Times-Bold",
    "Times-BoldItalic",
    "Helvetica",
    "Arial",
    "Helvetica-Oblique",
    "Helvetica-Bold",
    "Helvetica-BoldOblique",
    "Courier",
    "Courier-Oblique",
    "Courier-Bold",
    "Courier-BoldOblique",
    "Symbol",
    "Bookman-Light",
    "Bookman-Demi",
    "NewCenturySchlbk-Roman",
    "AvantGarde-Book",
    "Palatino-Roman",
    "ZapfChancery-MediumItalic",
    "ZapfDingbats",
]

sys_fonts = list_system_fonts()

println("\n==== GR 官方支持字体 ====")
foreach(f -> println("  ", f), gr_fonts)
println("\n==== 系统可用字体（部分） ====")
foreach(f -> println("  ", f), sys_fonts[1:min(end, 20)])  # 只显示前20个

# 输出交集和差集
common_fonts = intersect(gr_fonts, sys_fonts)
missing_fonts = setdiff(gr_fonts, sys_fonts)
println("\n==== 同时被 GR 和系统支持的字体 ====")
foreach(f -> println("  ", f), common_fonts)
println("\n==== GR 支持但系统未检测到的字体 ====")
foreach(f -> println("  ", f), missing_fonts)

# === 自动为每种可用字体生成测试图 ===
using Random
example_text = "The quick brown fox jumps over the lazy dog. 1234567890"
output_paths = String[]
for fname in common_fonts
    try
        x = 1:10
        y = rand(10)
        plt = plot(
            [1, 2],
            [1, 2],
            legend = false,
            framestyle = :none,
            grid = false,
            axis = false,
            title = fname,
            titlefont = font(32, fname),
            xlabel = "",
            ylabel = "",
        )
        annotate!(1.5, 1.5, Plots.text(example_text, 20, fname))
        outpath = joinpath(
            "scripts",
            "font_test_" * replace(fname, r"[^A-Za-z0-9]" => "_") * ".png",
        )
        savefig(plt, outpath)
        push!(output_paths, outpath)
    catch e
        @warn "生成字体测试图失败: $fname ($e)"
    end
end
println("\n已生成以下字体测试图：")
foreach(println, output_paths)
