import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/量子态与叠加原理_双缝干涉.py）

/// 双缝干涉强度分布的纯函数核。
enum DoubleSlitMath {

    /// 非相对论电子德布罗意波长 λ = h/√(2·m_e·E)。
    static func deBroglieWavelength(energyJoule E: Double, h: Double, me: Double) -> Double {
        h / (2 * me * E).squareRoot()
    }

    /// 干涉强度 I(x) = 4·I₀·cos²(π·d·x/(λL))。
    static func intensity(x: Double, lambda lam: Double, d: Double, L: Double,
                          I0: Double = 1) -> Double {
        let phase = .pi * d * x / (lam * L)
        return 4 * I0 * cos(phase) * cos(phase)
    }

    /// 条纹间距 Δx = λL/d。
    static func fringeSpacing(lambda lam: Double, d: Double, L: Double) -> Double {
        lam * L / d
    }
}

// MARK: - 模块

/// 笔记 09《量子态与叠加原理》：双缝干涉强度分布 + 条纹间距验证。
/// 实时档：4000 点解析曲线，电子能量/缝距/屏距拖动即重算。
struct DoubleSlitModule: SimModule {

    let meta = ModuleMeta(
        id: "量子态与叠加原理_双缝干涉", title: "量子态叠加 · 双缝干涉",
        subtitle: "I(x) = 4I₀cos²(πdx/λL)——叠加态的干涉条纹",
        category: .quantumStates, noteNumber: 9, tier: .realtime, difficulty: .basic,
        keywords: ["双缝", "干涉", "叠加", "superposition", "条纹间距",
                   "德布罗意", "de Broglie", "概率幅", "杨氏双缝"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "E", title: "电子动能", symbol: "E", unit: "eV",
                               range: 10...1000, defaultValue: 100,
                               scale: .log, decimalPlaces: 1)),
            .slider(SliderSpec(key: "d", title: "缝距", symbol: "d", unit: "nm",
                               range: 50...500, defaultValue: 100,
                               scale: .linear, decimalPlaces: 0)),
            .slider(SliderSpec(key: "L", title: "屏距", symbol: "L", unit: "m",
                               range: 0.1...2.0, defaultValue: 1.0,
                               scale: .linear, decimalPlaces: 2)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "双缝干涉强度 I(x) = 4I₀cos²(πdx/λL)",
                xAxis: .init(label: "屏位置 x (mm)"),
                yAxis: .init(label: "相对强度 I/I₀"),
                seriesNames: ["I(x)"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let eEv = input.slider("E")
        let d = input.slider("d") * 1e-9
        let L = input.slider("L")
        let h = try constants.value("h")
        let me = try constants.value("m_e")
        let e = try constants.value("e")

        let lam = DoubleSlitMath.deBroglieWavelength(energyJoule: eEv * e, h: h, me: me)
        let fringe = DoubleSlitMath.fringeSpacing(lambda: lam, d: d, L: L)

        // 与脚本一致：x ∈ [-0.02, 0.02] m × 4000 点
        let x = Num.linspace(-0.02, 0.02, count: 4000)
        let intensity = x.map {
            DoubleSlitMath.intensity(x: $0, lambda: lam, d: d, L: L)
        }

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: LineSeriesSpec(
                        xAxis: .init(label: "屏位置 x (mm)"),
                        yAxis: .init(label: "相对强度 I/I₀"),
                        seriesNames: ["I(x)"]),
                    series: [.init(name: "I(x)",
                                   points: Num.strided(x.map { $0 * 1e3 }, intensity, stride: 8))],
                    referenceLines: [
                        ReferenceLine(id: "f1", label: String(format: "+Δx = %.3f mm", fringe * 1e3),
                                      axis: .x, value: fringe * 1e3),
                        ReferenceLine(id: "fm1", label: String(format: "−Δx = −%.3f mm", fringe * 1e3),
                                      axis: .x, value: -fringe * 1e3),
                    ])),
            ],
            summary: [
                .init(id: "lam", title: "德布罗意波长 λ",
                      value: String(format: "%.4f Å", lam * 1e10), note: "h/√(2m_e·E)"),
                .init(id: "fringe", title: "条纹间距 Δx",
                      value: String(format: "%.4f mm", fringe * 1e3), note: "λL/d"),
                .init(id: "imax", title: "最大强度",
                      value: "4 I₀", note: "x = 0 中央主极大"),
                .init(id: "imin", title: "最小强度",
                      value: "0", note: "半整数条纹处"),
            ],
            theory: TheoryCard(
                title: "双缝干涉与叠加原理",
                formulas: [
                    "I(x) = 4I₀·cos²(π·d·x/(λL))（两缝概率幅相加后取模方）",
                    "条纹间距 Δx = λL/d（相邻极大间距）",
                    "λ = h/p = h/√(2m_e·E)（非相对论德布罗意波长）",
                    "叠加原理：|ψ₁+ψ₂|² = |ψ₁|² + |ψ₂|² + 2Re(ψ₁*ψ₂)",
                ],
                reading: "拖 E 看条纹变密（波长短）；拖 d 看条纹间距反向变化——Δx = λL/d 一目了然。"))
    }
}
