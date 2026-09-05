import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/物质波包_高斯波包演化.py）

/// 自由电子高斯波包解析演化的纯函数核。
enum WavePacketMath {

    /// 色散展宽宽度 σ(t) = σ₀·√(1 + (ℏt/(2mσ₀²))²)。
    static func sigmaT(_ t: Double, sigma0: Double, hbar: Double, m: Double) -> Double {
        sigma0 * (1 + pow(hbar * t / (2 * m * sigma0 * sigma0), 2)).squareRoot()
    }

    /// 归一化高斯概率密度 |ψ|² = exp(−(x−xc)²/2σ²) / √(2πσ²)。
    static func psi2(x: Double, xc: Double, sigma: Double) -> Double {
        exp(-pow(x - xc, 2) / (2 * sigma * sigma)) / (2 * .pi * sigma * sigma).squareRoot()
    }

    /// 群速度 v_g = ℏk₀/m（= 粒子速度 v）。
    static func groupVelocity(k0: Double, hbar: Double, m: Double) -> Double {
        hbar * k0 / m
    }

    /// 相速度 v_p = ℏk₀/(2m)（= v/2，仅相位、不载信息）。
    static func phaseVelocity(k0: Double, hbar: Double, m: Double) -> Double {
        hbar * k0 / (2 * m)
    }

    /// 数值宽度：√(∫ψ²(x−x̂)²dx / ∫ψ²dx)（均匀网格矩形求和，与脚本一致）。
    static func numericWidth(xs: [Double], psi2: [Double], xpeak: Double) -> Double {
        let dx = xs[1] - xs[0]
        var num = 0.0, den = 0.0
        for i in xs.indices {
            num += psi2[i] * pow(xs[i] - xpeak, 2) * dx
            den += psi2[i] * dx
        }
        return (num / den).squareRoot()
    }
}

// MARK: - 模块

/// 笔记 07《物质波包与波粒统一》：实验室系平移（ps）+ 共动系展宽（ns）双图组。
/// 实时档：解析曲线族，速度与初始宽度拖动即重算。
struct WavePacketModule: SimModule {

    /// 实验室系时刻（ps）与共动系时刻（ns），与脚本一致。
    static let labTimesPs = [0.0, 0.5, 1.0, 2.0]
    static let comTimesNs = [0.0, 0.1, 0.3, 0.5]

    let meta = ModuleMeta(
        id: "物质波包_高斯波包演化", title: "物质波包 · 高斯波包演化",
        subtitle: "群速度平移 vs 色散展宽——波粒统一的直观图像",
        category: .quantumStates, noteNumber: 7, tier: .realtime, difficulty: .basic,
        keywords: ["波包", "高斯", "群速度", "相速度", "色散", "展宽",
                   "德布罗意", "wave packet", "Gaussian"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "beta", title: "速度（c 的倍数）", symbol: "v/c",
                               unit: "c", range: 0.001...0.1, defaultValue: 0.01,
                               scale: .log, decimalPlaces: 3)),
            .slider(SliderSpec(key: "sigma0", title: "初始宽度", symbol: "σ₀", unit: "μm",
                               range: 0.05...1.0, defaultValue: 0.2,
                               scale: .linear, decimalPlaces: 2)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "实验室系：波包峰以 v_g = v 平移",
                xAxis: .init(label: "x (μm)"),
                yAxis: .init(label: "|ψ|²（归一）"),
                seriesNames: Self.labTimesPs.map { String(format: "t = %.1f ps", $0) })),
            .lineSeries(LineSeriesSpec(
                title: "共动系 x′ = x − v_g t：宽度按 σ(t) 展宽",
                xAxis: .init(label: "x′ (μm)"),
                yAxis: .init(label: "|ψ|²（归一）"),
                seriesNames: Self.comTimesNs.map { String(format: "t = %.1f ns", $0) })),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let beta = input.slider("beta")
        let sigma0 = input.slider("sigma0") * 1e-6
        let hbar = try constants.value("hbar")
        let me = try constants.value("m_e")
        let c = try constants.value("c")

        let v = beta * c
        let k0 = me * v / hbar
        let lam = 2 * .pi / k0
        let vg = WavePacketMath.groupVelocity(k0: k0, hbar: hbar, m: me)
        let vp = WavePacketMath.phaseVelocity(k0: k0, hbar: hbar, m: me)

        // 实验室系（ps 尺度，窗口 [-1,9] μm，2000 点与脚本一致）
        let xlab = Num.linspace(-1e-6, 9e-6, count: 2000)
        var labSeries: [SeriesPoints] = []
        for tps in Self.labTimesPs {
            let t = tps * 1e-12
            let xc = vg * t
            let st = WavePacketMath.sigmaT(t, sigma0: sigma0, hbar: hbar, m: me)
            let raw = xlab.map { WavePacketMath.psi2(x: $0, xc: xc, sigma: st) }
            let mx = raw.max() ?? 1
            labSeries.append(.init(
                name: String(format: "t = %.1f ps", tps),
                points: Num.strided(xlab.map { $0 * 1e6 }, raw.map { $0 / mx }, stride: 4)))
        }

        // 共动系（ns 尺度，窗口 [-1,1.5] μm）
        let xp = Num.linspace(-1e-6, 1.5e-6, count: 2000)
        var comSeries: [SeriesPoints] = []
        for tns in Self.comTimesNs {
            let t = tns * 1e-9
            let st = WavePacketMath.sigmaT(t, sigma0: sigma0, hbar: hbar, m: me)
            let raw = xp.map { WavePacketMath.psi2(x: $0, xc: 0, sigma: st) }
            let mx = raw.max() ?? 1
            comSeries.append(.init(
                name: String(format: "t = %.1f ns", tns),
                points: Num.strided(xp.map { $0 * 1e6 }, raw.map { $0 / mx }, stride: 4)))
        }

        let sigmaHalfNs = WavePacketMath.sigmaT(0.5e-9, sigma0: sigma0, hbar: hbar, m: me)

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: LineSeriesSpec(
                        xAxis: .init(label: "x (μm)"), yAxis: .init(label: "|ψ|²（归一）"),
                        seriesNames: Self.labTimesPs.map { String(format: "t = %.1f ps", $0) }),
                    series: labSeries)),
                .lineSeries(LineSeriesData(
                    spec: LineSeriesSpec(
                        xAxis: .init(label: "x′ (μm)"), yAxis: .init(label: "|ψ|²（归一）"),
                        seriesNames: Self.comTimesNs.map { String(format: "t = %.1f ns", $0) }),
                    series: comSeries)),
            ],
            summary: [
                .init(id: "lam", title: "德布罗意波长",
                      value: String(format: "%.2f pm", lam * 1e12), note: "2π/k₀"),
                .init(id: "vg", title: "群速度 v_g",
                      value: String(format: "%.4f c", vg / c), note: "= v（载信息）"),
                .init(id: "vp", title: "相速度 v_p",
                      value: String(format: "%.4f c", vp / c), note: "= v/2（仅相位）"),
                .init(id: "sig", title: "σ(0.5 ns)",
                      value: String(format: "%.3f μm", sigmaHalfNs * 1e6),
                      note: "σ₀√(1+(ℏt/2mσ₀²)²)"),
            ],
            theory: TheoryCard(
                title: "波包的两个速度",
                formulas: [
                    "σ(t) = σ₀·√(1 + (ℏt/(2mσ₀²))²)（色散展宽）",
                    "v_g = ℏk₀/m = v（波包整体以粒子速度平移）",
                    "v_p = ℏk₀/(2m) = v/2（相速度不载信息）",
                    "λ = 2π/k₀ = h/(m_e·v)（德布罗意波长）",
                ],
                reading: "左图看「粒子」——峰以 v 平移；右图看「波」——宽度随时间散开。调小 σ₀ 展宽更快（不确定性关系）。"))
    }
}
