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
            .slider(SliderSpec(key: "a", title: "结区边长", symbol: "a", unit: "µm",
                               range: 1...50, defaultValue: 10,
                               decimalPlaces: 1)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "BCS 序参量（能隙）随绝对温度",
                xAxis: .init(label: "T (K)", scale: .log),
                yAxis: .init(label: "Δ(T)/Δ(0)"),
                seriesNames: ["Δ(T)/Δ(0)"])),
            .lineSeries(LineSeriesSpec(
                title: "约瑟夫森结临界电流（Fraunhofer，绝对磁通）",
                xAxis: .init(label: "磁通密度 B (mT)"),
                yAxis: .init(label: "I_c(Φ)/I_c0"),
                seriesNames: ["I_c/I_c0"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let tc = input.slider("Tc")
        let aUm = input.slider("a")
        let kB = try constants.value("kB")
        let eV = try constants.value("eV")

        // W13 修正 B1：约化坐标 (T/T_c) 下曲线族恒为普适形状，Tc 扰动视觉 0%——
        // 改绝对温度 log 域 0.1…100 K（覆盖滑杆全域），能隙闭合点随 Tc 真实平移。
        let temps = Num.logspace(0.1, 100, count: 600)
        let ratio = temps.map { SuperconductivityMath.gapRatio($0 / tc) }

        // W13：Fraunhofer 同理改绝对磁通密度 B（结面积 A = a²，Φ = B·A）——
        // 场周期 Φ₀/A 随结区边长真实伸缩（结越大周期越小）。
        let phi0 = try constants.value("Phi0")  // 库珀对磁通量子
        let area = pow(aUm * 1e-6, 2)  // m²
        let phiRatio = Num.linspace(-3, 3, count: 900)
        let bGrid = phiRatio.map { $0 * phi0 / area * 1e3 }  // mT
        let ic = phiRatio.map(SuperconductivityMath.fraunhofer)

        let delta0 = 1.76 * kB * tc  // BCS: Δ(0) = 1.76 k_B T_c
        let ratioHalf = SuperconductivityMath.gapRatio(0.5)
        let periodMT = phi0 / area * 1e3  // 场周期（mT）

        let chart0 = LineSeriesData(
            spec: .init(xAxis: .init(label: "T (K)", scale: .log),
                        yAxis: .init(label: "Δ(T)/Δ(0)"),
                        seriesNames: ["Δ(T)/Δ(0)"]),
            series: [.init(name: "Δ(T)/Δ(0)",
                           points: Num.strided(temps, ratio, stride: 2))],
            referenceLines: [ReferenceLine(label: String(format: "T_c = %.2f K", tc),
                                           axis: .x, value: tc)],
            pointMarkers: [
                PointMarker(x: tc / 2, y: ratioHalf,
                            label: String(format: "Δ(0.5T_c)/Δ(0) = %.4f", ratioHalf)),
            ])

        let chart1 = LineSeriesData(
            spec: .init(xAxis: .init(label: "磁通密度 B (mT)"),
                        yAxis: .init(label: "I_c(Φ)/I_c0"),
                        seriesNames: ["I_c/I_c0"]),
            series: [.init(name: "I_c/I_c0",
                           points: Num.strided(bGrid, ic, stride: 2))],
            referenceLines: [
                ReferenceLine(label: "首个零点 Φ = Φ₀", axis: .x,
                              value: periodMT, style: .subtle),
                ReferenceLine(label: "−Φ₀", axis: .x, value: -periodMT, style: .subtle),
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
                .init(id: "period", title: "场周期 Φ₀/A",
                      value: String(format: "%.4f mT", periodMT),
                      note: String(format: "结区 %.1f×%.1f µm²（a 越大周期越小）", aUm, aUm)),
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
