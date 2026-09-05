import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/角动量与自旋_两自旋单态三重态.py）

/// Bell 态张量积、部分迹与纠缠熵的纯函数核（复数 re/im 分离小矩阵）。
enum TwoSpinMath {

    /// 内部小复数（避免为 4×4 矩阵引第三方复数库）。
    struct Z {
        var re: Double, im: Double
        static let zero = Z(re: 0, im: 0)
        static func *(lhs: Double, rhs: Z) -> Z { Z(re: lhs * rhs.re, im: lhs * rhs.im) }
        var abs2: Double { re * re + im * im }
    }

    /// Bell 态（基序 |↑↑⟩,|↑↓⟩,|↓↑⟩,|↓↓⟩，全部实系数）。
    /// |S⟩ = (|↑↓⟩−|↓↑⟩)/√2，|T+⟩ = |↑↑⟩，|T0⟩ = (|↑↓⟩+|↓↑⟩)/√2，|T−⟩ = |↓↓⟩。
    static func bellState(_ id: String) -> [Z] {
        let s2 = 2.0.squareRoot()
        switch id {
        case "Tp": return [Z(re: 1, im: 0), .zero, .zero, .zero]
        case "T0": return [.zero, Z(re: 1 / s2, im: 0), Z(re: 1 / s2, im: 0), .zero]
        case "Tm": return [.zero, .zero, .zero, Z(re: 1, im: 0)]
        default:   return [.zero, Z(re: 1 / s2, im: 0), Z(re: -1 / s2, im: 0), .zero]
        }
    }

    /// 密度矩阵 ρ = |ψ⟩⟨ψ|（4×4，行主序）。
    static func outer(_ psi: [Z]) -> [[Z]] {
        var rho = [[Z]](repeating: [Z](repeating: .zero, count: 4), count: 4)
        for i in 0..<4 {
            for j in 0..<4 {
                // ⟨ψ| 的分量是 conj(ψ)：ρ[i][j] = ψ[i]·conj(ψ[j])
                rho[i][j] = Z(re: psi[i].re * psi[j].re + psi[i].im * psi[j].im,
                              im: psi[i].im * psi[j].re - psi[i].re * psi[j].im)
            }
        }
        return rho
    }

    /// 对子系统 B（后一比特）部分迹 → 2×2 ρA：ρA[a,a'] = ρ[2a,2a'] + ρ[2a+1,2a'+1]。
    static func partialTraceOverB(_ rho: [[Z]]) -> [[Z]] {
        var r = [[Z]](repeating: [Z](repeating: .zero, count: 2), count: 2)
        for a in 0..<2 {
            for a2 in 0..<2 {
                let s0 = add(rho[2 * a][2 * a2], rho[2 * a + 1][2 * a2 + 1])
                r[a][a2] = s0
            }
        }
        return r
    }

    private static func add(_ a: Z, _ b: Z) -> Z { Z(re: a.re + b.re, im: a.im + b.im) }

    /// 2×2 厄米矩阵本征值（闭式，与 np.linalg.eigvalsh 对齐后排序）：
    /// λ = t/2 ± √((d/2)² + |m₁₂|²)，t = 迹，d = 对角差。
    static func hermitianEigenvalues(_ m: [[Z]]) -> [Double] {
        let t = m[0][0].re + m[1][1].re
        let d = (m[0][0].re - m[1][1].re) / 2
        let off = m[0][1].abs2.squareRoot()
        let r = (d * d + off * off).squareRoot()
        return [t / 2 - r, t / 2 + r].sorted { $0 < $1 }
    }

    /// von Neumann 熵（nat）：−Σ λ ln λ（过滤 λ ≤ 1e-12，脚本口径）。
    static func vonNeumannEntropy(_ m: [[Z]]) -> Double {
        var s = 0.0
        for ev in hermitianEigenvalues(m) where ev > 1e-12 {
            s -= ev * log(ev)
        }
        return s
    }

    /// 两自旋关联：单态 ⟨σ₁·â σ₂·b̂⟩ = −cos θ；三重态 T0 = +cos θ。
    static func correlation(theta: Double, singlet: Bool) -> Double {
        (singlet ? -1 : 1) * cos(theta)
    }
}

// MARK: - 模块

/// 笔记 16《角动量与自旋》§2.6/§6.1：两自旋 1/2 的 Bell 态——
/// 单态最大纠缠（ρA = 𝟙/2，S = ln 2），三重态 T0 同样纠缠，T± 为直积态；
/// 自旋关联反平行（−cos θ）vs 平行（+cos θ）。
struct TwoSpinModule: SimModule {

    let meta = ModuleMeta(
        id: "角动量与自旋_两自旋单态三重态", title: "两自旋 · 单态与三重态",
        subtitle: "张量积、部分迹与纠缠熵——Bell 基下的四个两比特态",
        category: .angularCentral, noteNumber: 16, tier: .seconds, difficulty: .basic,
        keywords: ["自旋", "单态", "三重态", "Bell 态", "张量积", "部分迹",
                   "密度矩阵", "纠缠熵", "spin", "singlet", "triplet"])

    var params: [ParamSpec] {
        [
            .discrete(DiscreteSpec(
                key: "state", title: "Bell 态",
                options: [
                    .init(id: "S", title: "单态 |S⟩", subtitle: "(|↑↓⟩−|↓↑⟩)/√2 · 最大纠缠"),
                    .init(id: "Tp", title: "三重态 |T+⟩", subtitle: "|↑↑⟩ · 直积态"),
                    .init(id: "T0", title: "三重态 |T0⟩", subtitle: "(|↑↓⟩+|↓↑⟩)/√2 · 最大纠缠"),
                    .init(id: "Tm", title: "三重态 |T−⟩", subtitle: "|↓↓⟩ · 直积态"),
                ],
                defaultOptionID: "S")),
        ]
    }

    var charts: [ChartSpec] {
        [
            .heatmap(HeatmapSpec(
                title: "约化密度矩阵 Re ρA",
                xAxis: .init(label: "A 基（列）"), yAxis: .init(label: "A 基（行）"),
                valueLabel: "Re ρ", diverging: true)),
            .lineSeries(LineSeriesSpec(
                title: "自旋关联 ⟨σ₁·â σ₂·b̂⟩",
                xAxis: .init(label: "夹角 θ (rad)"), yAxis: .init(label: "关联"),
                seriesNames: ["单态 |S⟩", "三重态 |T0⟩"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let stateID = input.discrete("state")
        let psi = TwoSpinMath.bellState(stateID)
        let rho = TwoSpinMath.outer(psi)
        let rhoA = TwoSpinMath.partialTraceOverB(rho)
        let entropy = TwoSpinMath.vonNeumannEntropy(rhoA)
        let trace = rhoA[0][0].re + rhoA[1][1].re
        let zBias = rhoA[0][0].re - rhoA[1][1].re

        // 关联曲线（脚本原样双线：单态与 T0，θ ∈ [0, π] 200 点）
        let theta = Num.linspace(0, .pi, count: 200)
        let singletPts = theta.map { Point(x: $0, y: TwoSpinMath.correlation(theta: $0, singlet: true)) }
        let tripletPts = theta.map { Point(x: $0, y: TwoSpinMath.correlation(theta: $0, singlet: false)) }

        let entangled = abs(entropy - log(2.0)) < 1e-6

        return SimResult(
            charts: [
                .heatmap(HeatmapData(
                    spec: HeatmapSpec(
                        title: "约化密度矩阵 Re ρA",
                        xAxis: .init(label: "A 基（列）"), yAxis: .init(label: "A 基（行）"),
                        valueLabel: "Re ρ", diverging: true),
                    xTicks: ["|↑⟩", "|↓⟩"], yTicks: ["|↑⟩", "|↓⟩"],
                    values: [[rhoA[0][0].re, rhoA[0][1].re], [rhoA[1][0].re, rhoA[1][1].re]])),
                .lineSeries(LineSeriesData(
                    spec: LineSeriesSpec(
                        title: "自旋关联 ⟨σ₁·â σ₂·b̂⟩",
                        xAxis: .init(label: "夹角 θ (rad)"), yAxis: .init(label: "关联"),
                        seriesNames: ["单态 |S⟩", "三重态 |T0⟩"]),
                    series: [
                        .init(name: "单态 |S⟩", points: singletPts),
                        .init(name: "三重态 |T0⟩", points: tripletPts, colorIndex: 1),
                    ])),
            ],
            summary: [
                .init(id: "ent", title: "纠缠熵 S(ρA)",
                      value: String(format: "%.4f", entropy), note: "nat"),
                .init(id: "ebit", title: "纠缠量",
                      value: String(format: "%.4f", entropy / log(2.0)), note: "ebit"),
                .init(id: "trace", title: "Tr ρA",
                      value: String(format: "%.4f", trace), note: "归一性"),
                .init(id: "type", title: "态类型",
                      value: entangled ? "最大纠缠（ρA = 𝟙/2）" : "直积态（ρA 纯态）",
                      note: "⟨σz⟩A = " + String(format: "%.1f", zBias)),
            ],
            theory: TheoryCard(
                title: "单态 · 三重态与纠缠",
                formulas: [
                    "|S⟩ = (|↑↓⟩−|↓↑⟩)/√2，|T±⟩ = |↑↑⟩/|↓↓⟩，|T0⟩ = (|↑↓⟩+|↓↑⟩)/√2",
                    "ρA = TrB(|ψ⟩⟨ψ|)：对 B 部分迹——每个自由度只看自己的一半",
                    "最大纠缠 ⇔ ρA = 𝟙/2 ⇔ S(ρA) = ln 2（1 ebit）",
                    "关联：单态 ⟨σ₁·â σ₂·b̂⟩ = −â·b̂（反平行），T0 为 +â·b̂（平行）",
                ],
                reading: "切到 |T±⟩ 看 ρA 塌成单点投影（熵 0），再回 |S⟩/|T0⟩ 看 𝟙/2 热图（熵 ln 2）——"
                    + "同样的两个自旋，是否纠缠完全由约化后的『混』决定。右图转 θ：单态关联永远反号，"
                    + "这正是 EPR 反关联的源头。"))
    }
}
