import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/康普顿散射_位移与反冲.py）

/// 康普顿位移与反冲电子动能的纯函数核。
enum ComptonMath {

    /// 康普顿波长 λ_C = h/(m_e·c)。
    static func comptonWavelength(h: Double, me: Double, c: Double) -> Double {
        h / (me * c)
    }

    /// 位移 Δλ = λ_C(1 − cosθ)。
    static func deltaLambda(theta: Double, lambdaC: Double) -> Double {
        lambdaC * (1 - cos(theta))
    }

    /// 反冲电子动能 K_e = hc/λ₀ − hc/(λ₀+Δλ)（能量守恒：光子损失即电子所得）。
    static func recoilKinetic(theta: Double, lambda0: Double,
                              lambdaC: Double, hc: Double) -> Double {
        let dl = deltaLambda(theta: theta, lambdaC: lambdaC)
        return hc / lambda0 - hc / (lambda0 + dl)
    }
}

// MARK: - 模块

/// 笔记 05《康普顿散射》：Δλ = λ_C(1−cosθ) 与反冲电子动能双图。
/// 实时档：400 点 θ 网格两条解析曲线，入射波长拖动即重算。
struct ComptonModule: SimModule {

    let meta = ModuleMeta(
        id: "康普顿散射_位移与反冲", title: "康普顿散射 · 位移与反冲",
        subtitle: "Δλ = λ_C(1−cosθ)；反冲电子动能 K_e(θ)",
        category: .oldQuantum, noteNumber: 5, tier: .realtime, difficulty: .basic,
        keywords: ["康普顿", "Compton", "散射", "位移", "反冲", "X 射线",
                   "光子", "波长", "Mo Kα"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "lam0", title: "入射波长", symbol: "λ₀", unit: "pm",
                               range: 1...200, defaultValue: 71,
                               scale: .log, decimalPlaces: 1)),
            .constant(ConstantSpec(key: "lambda_C", title: "λ_C", note: "康普顿波长 h/(m_e·c)")),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "康普顿位移 Δλ = λ_C(1−cosθ)",
                xAxis: .init(label: "散射角 θ (°)"),
                yAxis: .init(label: "Δλ (pm)"),
                seriesNames: ["Δλ"])),
            .lineSeries(LineSeriesSpec(
                title: "反冲电子动能 K_e(θ)",
                xAxis: .init(label: "散射角 θ (°)"),
                yAxis: .init(label: "K_e (eV)"),
                seriesNames: ["K_e"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let lam0pm = input.slider("lam0")
        let lam0 = lam0pm * 1e-12
        let h = try constants.value("h")
        let me = try constants.value("m_e")
        let c = try constants.value("c")
        let e = try constants.value("e")
        let lamC = ComptonMath.comptonWavelength(h: h, me: me, c: c)
        let hc = h * c

        // 与 Python theta = linspace(0, π, 400) 一致
        let thetas = Num.linspace(0, .pi, count: 400)
        let degs = thetas.map { $0 * 180 / .pi }
        let dlamPm = thetas.map { ComptonMath.deltaLambda(theta: $0, lambdaC: lamC) * 1e12 }
        let keEv = thetas.map {
            ComptonMath.recoilKinetic(theta: $0, lambda0: lam0, lambdaC: lamC, hc: hc) / e
        }

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: LineSeriesSpec(
                        xAxis: .init(label: "散射角 θ (°)"),
                        yAxis: .init(label: "Δλ (pm)"),
                        seriesNames: ["Δλ"]),
                    series: [.init(name: "Δλ", points: Num.strided(degs, dlamPm, stride: 1))],
                    referenceLines: [ReferenceLine(
                        id: "2lc", label: String(format: "2λ_C = %.2f pm", 2 * lamC * 1e12),
                        axis: .y, value: 2 * lamC * 1e12)])),
                .lineSeries(LineSeriesData(
                    spec: LineSeriesSpec(
                        xAxis: .init(label: "散射角 θ (°)"),
                        yAxis: .init(label: "K_e (eV)"),
                        seriesNames: ["K_e"]),
                    series: [.init(name: "K_e", points: Num.strided(degs, keEv, stride: 1))])),
            ],
            summary: [
                .init(id: "lc", title: "康普顿波长 λ_C",
                      value: String(format: "%.6f pm", lamC * 1e12), note: "h/(m_e·c)"),
                .init(id: "dl180", title: "Δλ(180°)",
                      value: String(format: "%.4f pm", 2 * lamC * 1e12), note: "= 2λ_C"),
                .init(id: "ke180", title: "K_e(180°)",
                      value: String(format: "%.1f eV", keEv.last ?? 0),
                      note: "背散射最大反冲"),
            ],
            theory: TheoryCard(
                title: "康普顿散射",
                formulas: [
                    "Δλ = λ_C(1 − cosθ)，λ_C = h/(m_e·c) = 2.426 pm",
                    "K_e = hc/λ₀ − hc/(λ₀ + Δλ)（能量守恒）",
                    "θ = π 背散射位移最大 = 2λ_C，与 λ₀ 无关",
                    "位移只取决于角度——光子「粒子性」的直接证据",
                ],
                reading: "拖 λ₀ 看反冲能量标尺变化（位移曲线纹丝不动）；90° 位移恰为一个康普顿波长。"))
    }
}
