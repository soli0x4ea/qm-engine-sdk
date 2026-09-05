import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/绘景变换与密度矩阵_两能级Rabi.py）

/// 两能级 Rabi 振荡的纯函数核（ħ = 1）：H = (ω/2)σz + (Ω/2)σx。
/// 薛定谔绘景用 Pauli 代数 2×2 矩阵指数闭式，相互作用绘景用 RK4，
/// 两条 ⟨σz⟩ 曲线验证绘景等价。
enum RabiMath {

    /// 薛绘景演化后的 ⟨σz⟩（闭式矩阵指数；ψ₀ = |↑⟩）：
    /// U_S(t) = cos(Ω_R t/2)·I − i·sin(Ω_R t/2)·(cφ σx + sφ σz)。
    static func sigmaZSchrodinger(t: Double, omega: Double, bigOmega: Double) -> Double {
        let oR = (omega * omega + bigOmega * bigOmega).squareRoot()
        let c = cos(oR * t / 2), s = sin(oR * t / 2)
        let cphi = bigOmega / oR, sphi = omega / oR
        // U_S|↑⟩ = (cos − i·s·sφ, −i·s·cφ) ⇒ ⟨σz⟩ = cos² + s²sφ² − s²cφ²
        return c * c + s * s * sphi * sphi - s * s * cphi * cphi
    }

    /// Rabi 频率 Ω_R = √(ω² + Ω²)。
    static func rabiFrequency(omega: Double, bigOmega: Double) -> Double {
        (omega * omega + bigOmega * bigOmega).squareRoot()
    }

    /// 相互作用绘景 RHS：dψ/dt = −i·V_I(t)ψ，V_I = (Ω/2)(cos ωt·σx − sin ωt·σy)。
    /// ψ = (a, b) 复数；展开得（e^{∓iωt} 已内联）：
    /// ȧ = (Ω/2)(sin ωt·b_r + cos ωt·b_i) + i·(Ω/2)(sin ωt·b_i − cos ωt·b_r)，等。
    private static func rhs(
        t: Double, ar: Double, ai: Double, br: Double, bi: Double,
        omega: Double, bigOmega: Double
    ) -> (ar: Double, ai: Double, br: Double, bi: Double) {
        let s = sin(omega * t), c = cos(omega * t)
        let h = bigOmega / 2
        return (h * (s * br + c * bi),
                h * (s * bi - c * br),
                h * (-s * ar + c * ai),
                h * (-s * ai - c * ar))
    }

    /// RK4 一步（脚本逐句对应）。
    static func rk4Step(
        t: Double, _ psi: (ar: Double, ai: Double, br: Double, bi: Double),
        dt: Double, omega: Double, bigOmega: Double
    ) -> (ar: Double, ai: Double, br: Double, bi: Double) {
        let k1 = rhs(t: t, ar: psi.ar, ai: psi.ai, br: psi.br, bi: psi.bi,
                     omega: omega, bigOmega: bigOmega)
        let k2 = rhs(t: t + dt / 2,
                     ar: psi.ar + dt / 2 * k1.ar, ai: psi.ai + dt / 2 * k1.ai,
                     br: psi.br + dt / 2 * k1.br, bi: psi.bi + dt / 2 * k1.bi,
                     omega: omega, bigOmega: bigOmega)
        let k3 = rhs(t: t + dt / 2,
                     ar: psi.ar + dt / 2 * k2.ar, ai: psi.ai + dt / 2 * k2.ai,
                     br: psi.br + dt / 2 * k2.br, bi: psi.bi + dt / 2 * k2.bi,
                     omega: omega, bigOmega: bigOmega)
        let k4 = rhs(t: t + dt,
                     ar: psi.ar + dt * k3.ar, ai: psi.ai + dt * k3.ai,
                     br: psi.br + dt * k3.br, bi: psi.bi + dt * k3.bi,
                     omega: omega, bigOmega: bigOmega)
        return (psi.ar + dt / 6 * (k1.ar + 2 * k2.ar + 2 * k3.ar + k4.ar),
                psi.ai + dt / 6 * (k1.ai + 2 * k2.ai + 2 * k3.ai + k4.ai),
                psi.br + dt / 6 * (k1.br + 2 * k2.br + 2 * k3.br + k4.br),
                psi.bi + dt / 6 * (k1.bi + 2 * k2.bi + 2 * k3.bi + k4.bi))
    }

    /// 相互作用绘景全程扫描：t ∈ [0, 8π] 400 点，返回每点 ⟨σz⟩_I。
    static func sigmaZInteraction(omega: Double, bigOmega: Double)
        -> (ts: [Double], expZ: [Double], maxDiff: Double) {
        let ts = Num.linspace(0, 8 * .pi, count: 400)
        let dt = ts[1] - ts[0]
        var psi = (ar: 1.0, ai: 0.0, br: 0.0, bi: 0.0)   // ψ₀ = |↑⟩
        var expZ = [Double](repeating: 0, count: ts.count)
        expZ[0] = psi.ar * psi.ar + psi.ai * psi.ai - psi.br * psi.br - psi.bi * psi.bi
        var maxDiff = 0.0
        for n in 1..<ts.count {
            psi = rk4Step(t: ts[n - 1], psi, dt: dt, omega: omega, bigOmega: bigOmega)
            expZ[n] = psi.ar * psi.ar + psi.ai * psi.ai - psi.br * psi.br - psi.bi * psi.bi
            let diff = abs(expZ[n] - sigmaZSchrodinger(t: ts[n], omega: omega, bigOmega: bigOmega))
            if diff > maxDiff { maxDiff = diff }
        }
        return (ts, expZ, maxDiff)
    }
}

// MARK: - 模块

/// 笔记 18《绘景变换与密度矩阵》§2.4/§6.1：两能级 Rabi 振荡——
/// 薛定谔绘景（矩阵指数闭式）与相互作用绘景（RK4）给出完全相同的 ⟨σz⟩。
struct RabiModule: SimModule {

    let meta = ModuleMeta(
        id: "绘景变换与密度矩阵_两能级Rabi", title: "两能级 Rabi · 绘景等价",
        subtitle: "同一物理，两种坐标——薛定谔绘景闭式 vs 相互作用绘景 RK4",
        category: .formalTheory, noteNumber: 18, tier: .seconds, difficulty: .basic,
        keywords: ["Rabi 振荡", "绘景变换", "相互作用绘景", "RK4", "矩阵指数",
                   "两能级", "Rabi oscillation", "picture"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "omega", title: "失谐 / 进动率", symbol: "ω",
                               unit: "ħ", range: 0...3, defaultValue: 1, decimalPlaces: 2)),
            .slider(SliderSpec(key: "bigOmega", title: "驱动强度", symbol: "Ω",
                               unit: "ħ", range: 0.5...4, defaultValue: 2, decimalPlaces: 2)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "⟨σz⟩(t)：两种绘景的同一预言",
                xAxis: .init(label: "时间 t (ħ=1)"), yAxis: .init(label: "⟨σz⟩"),
                seriesNames: ["薛定谔绘景", "相互作用绘景 (RK4)"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let omega = input.slider("omega")
        let bigOmega = input.slider("bigOmega")

        let (ts, expI, maxDiff) = RabiMath.sigmaZInteraction(omega: omega, bigOmega: bigOmega)
        let expS = ts.map { RabiMath.sigmaZSchrodinger(t: $0, omega: omega, bigOmega: bigOmega) }

        let oR = RabiMath.rabiFrequency(omega: omega, bigOmega: bigOmega)
        // 翻转深度：⟨σz⟩ 极值 = (ω² − Ω²)/(ω² + Ω²)（共振时到 −1）
        let flipExtremum = (omega * omega - bigOmega * bigOmega)
            / (omega * omega + bigOmega * bigOmega)

        let sPts = ts.enumerated().map { Point(x: $0.element, y: expS[$0.offset]) }
        let iPts = ts.enumerated().map { Point(x: $0.element, y: expI[$0.offset]) }

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: LineSeriesSpec(
                        title: "⟨σz⟩(t)：两种绘景的同一预言",
                        xAxis: .init(label: "时间 t (ħ=1)"), yAxis: .init(label: "⟨σz⟩"),
                        seriesNames: ["薛定谔绘景", "相互作用绘景 (RK4)"]),
                    series: [
                        .init(name: "薛定谔绘景", points: sPts),
                        .init(name: "相互作用绘景 (RK4)", points: iPts, colorIndex: 1),
                    ])),
            ],
            summary: [
                .init(id: "oR", title: "Rabi 频率 Ω_R",
                      value: String(format: "%.4f", oR), note: "√(ω²+Ω²)"),
                .init(id: "period", title: "振荡周期",
                      value: String(format: "%.4f", 2 * .pi / oR), note: "2π/Ω_R"),
                .init(id: "flip", title: "翻转深度 ⟨σz⟩min",
                      value: String(format: "%.4f", flipExtremum),
                      note: omega < 1e-9 ? "共振完全翻转" : "失谐 ω ≠ 0 抑制翻转"),
                .init(id: "diff", title: "绘景偏差 max|Δ|",
                      value: String(format: "%.2e", maxDiff), note: "RK4 步长误差"),
            ],
            theory: TheoryCard(
                title: "绘景变换：换坐标不换物理",
                formulas: [
                    "H = (ω/2)σz + (Ω/2)σx，U_S(t) = exp(−iHt)（2×2 矩阵指数闭式）",
                    "相互作用绘景：ψ_I = e^{iH₀t}ψ_S，dψ_I/dt = −iV_I(t)ψ_I（RK4 数值积分）",
                    "⟨σz⟩(t) = cos²(Ω_R t/2) + sin²(Ω_R t/2)·(ω²−Ω²)/(ω²+Ω²)",
                    "两条曲线重合 ⇒ 绘景等价：可观测量与绘景选择无关",
                ],
                reading: "红实线（薛绘景闭式）与蓝虚线（相互作用绘景 RK4）应当完全重叠——"
                    + "绘景只是坐标变换，预言必须一致。拖大 ω 看失谐抑制翻转（振幅到不了 −1）；"
                    + "ω → 0 时共振完全翻转，周期 2π/Ω_R。"))
    }
}
