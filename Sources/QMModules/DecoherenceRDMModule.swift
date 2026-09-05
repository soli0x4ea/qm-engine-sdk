import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/绘景变换与密度矩阵_退相干约化密度矩阵.py）

/// Bell 单态经退相位噪声的约化密度矩阵纯函数核。
/// ρ_A(t) = ½[[1, e^(−Γt)], [e^(−Γt), 1]]——非对角元按 e^(−Γt) 衰减。
enum DecoherenceRDMMath {

    /// 非对角元（相干幅度）。
    static func offDiagonal(_ t: Double, gamma: Double) -> Double {
        exp(-gamma * t)
    }

    /// 纯度 Tr(ρ²) = (1 + e^(−2Γt))/2。
    static func purity(_ t: Double, gamma: Double) -> Double {
        0.5 * (1 + pow(offDiagonal(t, gamma: gamma), 2))
    }

    /// von Neumann 熵（ebit）：本征值 ½(1 ± e^(−Γt))，clip 1e-15 防端点 log(0)。
    static func entropyEbits(_ t: Double, gamma: Double) -> Double {
        let off = offDiagonal(t, gamma: gamma)
        let ev1 = max(0.5 * (1 + off), 1e-15)
        let ev2 = max(0.5 * (1 - off), 1e-15)
        let s = -(ev1 * log(ev1) + ev2 * log(ev2))
        return s / log(2)
    }
}

// MARK: - 模块

/// 笔记 18《绘景变换与密度矩阵》§2.6/§6.2：纠缠单态退相干——
/// 非对角元 e^(−Γt) 衰减、纯度跌落、熵升至 1 ebit。实时档。
struct DecoherenceRDMModule: SimModule {

    let meta = ModuleMeta(
        id: "绘景变换与密度矩阵_退相干约化密度矩阵",
        title: "退相干 · 约化密度矩阵",
        subtitle: "非对角元 e^(−Γt) 衰减 + 纯度与纠缠熵",
        category: .formalTheory, noteNumber: 18, tier: .realtime, difficulty: .basic,
        keywords: ["退相干", "dephasing", "密度矩阵", "约化", "纯度", "purity",
                   "熵", "entropy", "ebit", "纠缠", "Bell", "单态"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "Gamma", title: "退相位速率", symbol: "Γ", unit: "1/t",
                               range: 0.1...4, defaultValue: 1, decimalPlaces: 2)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "退相干：非对角元衰减与纯度跌落",
                xAxis: .init(label: "t"),
                yAxis: .init(label: "幅度"),
                seriesNames: ["ρ₁₂(t)", "Tr ρ²"])),
            .lineSeries(LineSeriesSpec(
                title: "纠缠熵升至 1 ebit",
                xAxis: .init(label: "t"),
                yAxis: .init(label: "熵 (ebit)"),
                seriesNames: ["S(ρ_A)/ln2"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let gamma = input.slider("Gamma")
        // 与 Python t = linspace(0, 4, 300) 一致
        let t = Num.linspace(0, 4, count: 300)

        let offdiag = t.map { DecoherenceRDMMath.offDiagonal($0, gamma: gamma) }
        let pur = t.map { DecoherenceRDMMath.purity($0, gamma: gamma) }
        let entropy = t.map { DecoherenceRDMMath.entropyEbits($0, gamma: gamma) }

        let chart0 = LineSeriesData(
            spec: .init(xAxis: .init(label: "t"), yAxis: .init(label: "幅度"),
                        seriesNames: ["ρ₁₂(t)", "Tr ρ²"]),
            series: [.init(name: "ρ₁₂(t)", points: Num.strided(t, offdiag, stride: 2)),
                     .init(name: "Tr ρ²", points: Num.strided(t, pur, stride: 2))])
        let chart1 = LineSeriesData(
            spec: .init(xAxis: .init(label: "t"), yAxis: .init(label: "熵 (ebit)"),
                        seriesNames: ["S(ρ_A)/ln2"]),
            series: [.init(name: "S(ρ_A)/ln2", points: Num.strided(t, entropy, stride: 2))],
            referenceLines: [ReferenceLine(label: "S = 1 ebit", axis: .y, value: 1.0,
                                          style: .subtle)])

        let tHalf = log(2) / gamma
        return SimResult(
            charts: [.lineSeries(chart0), .lineSeries(chart1)],
            summary: [
                .init(id: "gamma", title: "Γ", value: String(format: "%.2f 1/t", gamma),
                      note: "退相位速率"),
                .init(id: "halflife", title: "相干半衰期",
                      value: String(format: "%.2f t", tHalf),
                      note: "ln2/Γ"),
                .init(id: "asymptote", title: "t→∞",
                      value: "纯度 → 0.5，熵 → 1 ebit", note: "最大混合态"),
            ],
            theory: TheoryCard(
                title: "退相位下的约化密度矩阵",
                formulas: [
                    "ρ_A(t) = ½[[1, e^(−Γt)], [e^(−Γt), 1]]（Bell 单态 + 退相位噪声）",
                    "纯度 Tr ρ² = (1 + e^(−2Γt))/2：1 → ½",
                    "熵 S = −Tr ρ ln ρ（ebit）：0 → 1（最大混合）",
                ],
                reading: "非对角元（相干）最先消亡而布居不变——退相干 ≠ 弛豫；"
                    + "熵升至 1 ebit 即初态纠缠被环境完全稀释。"))
    }
}
