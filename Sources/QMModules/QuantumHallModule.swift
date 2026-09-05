import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/量子霍尔_Landau能级与平台.py）

/// 朗道能级 + IQHE 平台玩具模型的纯函数核。
enum QuantumHallMath {

    /// 朗道能级 E_n(B) = (n + ½)ℏeB/m*，J（自由电子 m* = m_e）。
    static func landauLevel(_ B: Double, n: Int, hbar: Double, me: Double,
                            e: Double) -> Double {
        (Double(n) + 0.5) * hbar * e * B / me
    }

    /// 填充因子 ν = n_e·h/(eB)。
    static func fillingFactor(_ B: Double, ne: Double, h: Double, e: Double) -> Double {
        ne * h / (e * B)
    }

    /// 填充 Landau 能级数 N = max(1, ⌊ν⌋)。
    static func filledLevels(_ nu: Double) -> Double {
        max(1, nu.rounded(.down))
    }

    /// 霍尔电阻率平台值 ρ_xy = h/(Ne²)，Ω。
    static func rhoXY(_ B: Double, ne: Double, h: Double, e: Double) -> Double {
        let nu = fillingFactor(B, ne: ne, h: h, e: e)
        return h / (filledLevels(nu) * e * e)
    }

    /// 纵向电阻率（玩具模型）：距最近整数填充的距离 ×1000，0.45 截断。
    static func rhoXX(_ B: Double, ne: Double, h: Double, e: Double) -> Double {
        let nu = fillingFactor(B, ne: ne, h: h, e: e)
        let dist = abs(nu - nu.rounded())
        return 1000 * min(max(dist / 0.45, 0), 1)
    }
}

// MARK: - 模块

/// 笔记 28《量子霍尔效应》：朗道能级谱 E_n(B)（ℏω_c 线性）
/// + 整数霍尔平台 ρ_xy 与 ρ_xx dip（双 Y 轴）。实时档。
struct QuantumHallModule: SimModule {

    let meta = ModuleMeta(
        id: "量子霍尔_Landau能级与平台", title: "量子霍尔 · Landau 能级与平台",
        subtitle: "E_n = (n+½)ℏω_c + ρ_xy 平台 / ρ_xx dip（双 Y 轴）",
        category: .manyBody, noteNumber: 28, tier: .realtime, difficulty: .basic,
        keywords: ["量子霍尔", "Landau", "朗道能级", "霍尔", "平台", "von Klitzing",
                   "克里青常数", "填充因子", "二维电子气"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "ne", title: "电子面密度", symbol: "n_e",
                               unit: "m⁻²", range: 5e14...5e16, defaultValue: 2.4e15,
                               scale: .log, decimalPlaces: 2)),
            .constant(ConstantSpec(key: "RK", title: "R_K", note: "克里青常数 h/e²")),
            .constant(ConstantSpec(key: "hbar", title: "ℏ", note: "约化普朗克常量")),
            .constant(ConstantSpec(key: "m_e", title: "m_e", note: "自由电子质量")),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "朗道能级谱随磁场",
                xAxis: .init(label: "磁场 B (T)"),
                yAxis: .init(label: "朗道能级 E_n (meV)"),
                seriesNames: ["n=0", "n=1", "n=2", "n=3", "n=4", "n=5"])),
            .dualAxisLineSeries(DualAxisLineSeriesSpec(
                title: "整数霍尔：ρ_xy 平台与 ρ_xx dip",
                xAxis: .init(label: "磁场 B (T)"),
                primaryAxis: .init(label: "ρ_xy/R_K"),
                secondaryAxis: .init(label: "ρ_xx (a.u.)"),
                primaryNames: ["ρ_xy/R_K"],
                secondaryNames: ["ρ_xx (dips)"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let ne = input.slider("ne")
        let hbar = try constants.value("hbar")
        let me = try constants.value("m_e")
        let e = try constants.value("e")
        let h = try constants.value("h")
        let eV = try constants.value("eV")
        let rk = h / (e * e)

        // ---- 段 1：朗道能级（B = linspace(0, 15, 400)，与 Python 一致）----
        let B = Num.linspace(0, 15, count: 400)
        var landauSeries: [SeriesPoints] = []
        for n in 0..<6 {
            let en = B.map {
                QuantumHallMath.landauLevel($0, n: n, hbar: hbar, me: me, e: e) / eV * 1e3
            }
            landauSeries.append(.init(name: "n=\(n)", points: zip(B, en).map {
                Point(x: $0, y: $1)
            }))
        }

        // ---- 段 2：IQHE 玩具模型（B2 = linspace(1, 10, 1200)）----
        let B2 = Num.linspace(1, 10, count: 1200)
        let rhoXY = B2.map { QuantumHallMath.rhoXY($0, ne: ne, h: h, e: e) / rk }
        let rhoXX = B2.map { QuantumHallMath.rhoXX($0, ne: ne, h: h, e: e) }

        let chart0 = LineSeriesData(
            spec: .init(xAxis: .init(label: "磁场 B (T)"),
                        yAxis: .init(label: "朗道能级 E_n (meV)"),
                        seriesNames: ["n=0", "n=1", "n=2", "n=3", "n=4", "n=5"]),
            series: landauSeries)

        let chart1 = DualAxisLineSeriesData(
            spec: .init(xAxis: .init(label: "磁场 B (T)"),
                        primaryAxis: .init(label: "ρ_xy/R_K"),
                        secondaryAxis: .init(label: "ρ_xx (a.u.)"),
                        primaryNames: ["ρ_xy/R_K"],
                        secondaryNames: ["ρ_xx (dips)"]),
            primary: [.init(name: "ρ_xy/R_K",
                            points: Num.strided(B2, rhoXY, stride: 2))],
            secondary: [.init(name: "ρ_xx (dips)",
                              points: Num.strided(B2, rhoXX, stride: 2))])

        let hbarWc10 = hbar * e * 10 / me / eV * 1e3
        let degen10 = e * 10 / h
        let rhoAt5 = QuantumHallMath.rhoXY(5.0, ne: ne, h: h, e: e)

        return SimResult(
            charts: [.lineSeries(chart0), .dualAxisLineSeries(chart1)],
            summary: [
                .init(id: "rk", title: "R_K = h/e²",
                      value: String(format: "%.1f Ω", rk),
                      note: "2019 SI 精确值 25812.80745"),
                .init(id: "wc", title: "ℏω_c @10 T",
                      value: String(format: "%.3f meV", hbarWc10),
                      note: "自由电子回旋能"),
                .init(id: "degen", title: "简并度 eB/h @10 T",
                      value: String(format: "%.3e m⁻²", degen10),
                      note: "每能级每单位面积态数"),
                .init(id: "plateau", title: "ρ_xy(5 T)",
                      value: String(format: "%.1f Ω", rhoAt5),
                      note: String(format: "= R_K/%.0f", rk / rhoAt5)),
            ],
            theory: TheoryCard(
                title: "Landau 量子化与霍尔平台",
                formulas: [
                    "E_n(B) = (n + ½)ℏω_c，ω_c = eB/m：能级随 B 线性上移",
                    "简并度（单位面积）= eB/h：B 增大能级容量同步增大",
                    "填充因子 ν = n_e·h/(eB)：费米能级随 B 增大逐级穿越朗道能级",
                    "平台：ν 为整数时费米能级在能隙内 → ρ_xx → 0，ρ_xy = R_K/ν",
                ],
                reading: "主图 6 条朗道线随 B 抬升；副图 ρ_xy 呈阶梯（1/ν 量化）"
                    + "而 ρ_xx 在平台中心归零——费米能级在能隙内的那段就是平台宽度。"))
    }
}
