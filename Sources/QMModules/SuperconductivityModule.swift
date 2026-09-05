import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/超导_BCS与约瑟夫森.py）

/// BCS 能隙插值 + 约瑟夫森 Fraunhofer 调制的纯函数核。
enum SuperconductivityMath {

    /// BCS 能隙插值：Δ(T)/Δ(0) = tanh(1.76·√(T_c/T − 1))（T < T_c，否则 0）。
    static func gapRatio(_ tRatio: Double) -> Double {
        guard tRatio < 1 else { return 0 }
        return tanh(1.76 * sqrt(1 / tRatio - 1))
    }

    /// 约瑟夫森临界电流 I_c/I_c0 = |sinc(Φ/Φ₀)| = |sin(πx)/(πx)|；x=0 取 1。
    static func fraunhofer(_ phiRatio: Double) -> Double {
        let x = .pi * phiRatio
        if abs(phiRatio) < 1e-12 { return 1 }
        return abs(sin(x) / x)
    }
}

// MARK: - 模块

/// 笔记 27《超导电性与宏观量子现象》：BCS 能隙随温度（插值公式）
/// + 约瑟夫森结临界电流的 Fraunhofer 调制。实时档。
struct SuperconductivityModule: SimModule {

    let meta = ModuleMeta(
        id: "超导_BCS与约瑟夫森", title: "超导 · BCS 与约瑟夫森",
        subtitle: "能隙插值 tanh(1.76√(T_c/T−1)) + Fraunhofer 调制",
        category: .manyBody, noteNumber: 27, tier: .realtime, difficulty: .basic,
        keywords: ["超导", "BCS", "库珀对", "能隙", "临界温度", "Tc", "约瑟夫森",
                   "Josephson", "Fraunhofer", "磁通量子", "sinc"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "Tc", title: "临界温度", symbol: "T_c", unit: "K",
                               range: 0.5...100, defaultValue: 9.25,
                               scale: .log, decimalPlaces: 2)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "BCS 序参量（能隙）随温度",
                xAxis: .init(label: "T/T_c"),
                yAxis: .init(label: "Δ(T)/Δ(0)"),
                seriesNames: ["Δ(T)/Δ(0)"])),
            .lineSeries(LineSeriesSpec(
                title: "约瑟夫森结临界电流（Fraunhofer）",
                xAxis: .init(label: "磁通 Φ/Φ₀"),
                yAxis: .init(label: "I_c(Φ)/I_c0"),
                seriesNames: ["I_c/I_c0"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let tc = input.slider("Tc")
        let kB = try constants.value("kB")
        let eV = try constants.value("eV")

        // 与 Python Tr = linspace(0.02, 1.4, 600) 一致
        let tRatio = Num.linspace(0.02, 1.4, count: 600)
        let ratio = tRatio.map(SuperconductivityMath.gapRatio)

        // 与 Python Phi_ratio = linspace(−3, 3, 900) 一致
        let phi = Num.linspace(-3, 3, count: 900)
        let ic = phi.map(SuperconductivityMath.fraunhofer)

        let delta0 = 1.76 * kB * tc  // BCS: Δ(0) = 1.76 k_B T_c
        let ratioHalf = SuperconductivityMath.gapRatio(0.5)

        let chart0 = LineSeriesData(
            spec: .init(xAxis: .init(label: "T/T_c"),
                        yAxis: .init(label: "Δ(T)/Δ(0)"),
                        seriesNames: ["Δ(T)/Δ(0)"]),
            series: [.init(name: "Δ(T)/Δ(0)",
                           points: Num.strided(tRatio, ratio, stride: 2))],
            referenceLines: [ReferenceLine(label: "T/T_c = 1", axis: .x, value: 1)])

        let chart1 = LineSeriesData(
            spec: .init(xAxis: .init(label: "磁通 Φ/Φ₀"),
                        yAxis: .init(label: "I_c(Φ)/I_c0"),
                        seriesNames: ["I_c/I_c0"]),
            series: [.init(name: "I_c/I_c0",
                           points: Num.strided(phi, ic, stride: 2))],
            referenceLines: [
                ReferenceLine(label: "首个零点 Φ/Φ₀ = 1", axis: .x, value: 1, style: .subtle),
                ReferenceLine(label: "−1", axis: .x, value: -1, style: .subtle),
            ])

        return SimResult(
            charts: [.lineSeries(chart0), .lineSeries(chart1)],
            summary: [
                .init(id: "d0", title: "Δ(0)",
                      value: String(format: "%.3f meV", delta0 / eV * 1e3),
                      note: String(format: "1.76 k_B T_c（T_c = %.2f K）", tc)),
                .init(id: "half", title: "Δ(0.5T_c)/Δ(0)",
                      value: String(format: "%.4f", ratioHalf),
                      note: "低温平台——库珀对束缚坚固"),
                .init(id: "fraun", title: "I_c(0)",
                      value: String(format: "%.4f I_c0", 1.0),
                      note: "Fraunhofer 中心峰；零点在整数 Φ₀"),
            ],
            theory: TheoryCard(
                title: "BCS 能隙与约瑟夫森衍射",
                formulas: [
                    "Δ(T)/Δ(0) = tanh(1.76√(T_c/T − 1))：T_c 处闭合",
                    "Δ(0) = 1.76 k_BT_c（弱耦合 BCS 极限）",
                    "I_c(Φ) = I_c0·|sin(πΦ/Φ₀)/(πΦ/Φ₀)|——结区磁通的衍射图样",
                    "Φ₀ = h/2e ≈ 2.07×10⁻¹⁵ Wb（库珀对磁通量子）",
                ],
                reading: "能隙在 T→0 近乎平台、在 T_c 附近按 √(1−T/T_c) 收口；"
                    + "Fraunhofer 图样把结当成「量子衍射光栅」——磁通每加一个 Φ₀，"
                    + "临界电流过零一次。"))
    }
}
