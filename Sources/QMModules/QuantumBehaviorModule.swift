import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/量子行为_三实验对照.py）

/// 量子行为四分布对照的纯函数核。自然单位 λ = L = 1（仅几何比例有意义）。
enum QuantumBehaviorMath {

    /// np.sinc：归一化 sinc（sin(πx)/(πx)，sinc(0)=1）。
    static func sinc(_ x: Double) -> Double {
        x == 0 ? 1 : sin(.pi * x) / (.pi * x)
    }

    /// 单缝远场幅的包络（夫琅禾费近似下不随缝位移）与相位。
    /// shift = ±d/2 → 相位 exp(∓iπ·d·x/(λL))（s = sign(shift)，shift=0 取 s=1）。
    static func slitAmplitude(x: Double, shift: Double,
                              lam: Double = 1, L: Double = 1,
                              a: Double, d: Double) -> (re: Double, im: Double) {
        let envelope = sinc(a * x / (lam * L))
        let s = shift != 0 ? (shift > 0 ? 1.0 : -1.0) : 1.0
        let phase = -s * .pi * d * x / (lam * L)
        return (envelope * cos(phase), envelope * sin(phase))
    }

    /// 四分布在同一 x 网格上求值（归一化逻辑与 Python main() 一致）：
    /// 返回 (x, bullets, waveB, electron, watched)。
    /// - bullets: P1+P2（未单独归一化，整体归一）
    /// - waveB:   P2 单独归一（两独立源之一，Python 画的是 P2n）
    /// - electron: |φ1+φ2|² 归一
    /// - watched: P1n + P2n（各自归一后相加）
    static func distributions(d: Double, a: Double,
                              count: Int = 3000
    ) -> (x: [Double], bullets: [Double], waveB: [Double],
          electron: [Double], watched: [Double]) {
        let x = Num.linspace(-3, 3, count: count)
        let dx = x[1] - x[0]
        var p1 = [Double](), p2 = [Double](), pe = [Double]()
        p1.reserveCapacity(count); p2.reserveCapacity(count); pe.reserveCapacity(count)
        for xi in x {
            let f1 = slitAmplitude(x: xi, shift: -d / 2, a: a, d: d)
            let f2 = slitAmplitude(x: xi, shift: +d / 2, a: a, d: d)
            let a1 = f1.re * f1.re + f1.im * f1.im
            let a2 = f2.re * f2.re + f2.im * f2.im
            p1.append(a1); p2.append(a2)
            let sr = f1.re + f2.re, si = f1.im + f2.im
            pe.append(sr * sr + si * si)
        }
        let s1 = p1.reduce(0, +) * dx, s2 = p2.reduce(0, +) * dx
        let se = pe.reduce(0, +) * dx
        let bullets = zip(p1, p2).map(+)
        let sb = bullets.reduce(0, +) * dx
        return (x,
                bullets.map { $0 / sb },
                p2.map { $0 / s2 },
                pe.map { $0 / se },
                zip(p1, p2).map { $0.0 / s1 + $0.1 / s2 })
    }
}

// MARK: - 模块

/// 笔记 01《量子行为》：子弹 / 水波 / 电子 / 被监视电子 四分布对照（费曼第三卷第 1 章）。
/// 实时档：3000 点网格四条曲线，缝距/缝宽拖动即重算。
struct QuantumBehaviorModule: SimModule {

    let meta = ModuleMeta(
        id: "量子行为_三实验对照", title: "量子行为 · 三实验对照",
        subtitle: "子弹 / 水波 / 电子 / 被监视电子——同几何三种行为",
        category: .oldQuantum, noteNumber: 1, tier: .realtime, difficulty: .basic,
        keywords: ["双缝", "干涉", "费曼", "Feynman", "概率幅", "which-way",
                   "夫琅禾费", "子弹", "水波", "电子", "监视"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "d", title: "缝距", symbol: "d", unit: "λL",
                               range: 0.2...3.0, defaultValue: 1.0,
                               scale: .linear, decimalPlaces: 2)),
            .slider(SliderSpec(key: "a", title: "缝宽", symbol: "a", unit: "λL",
                               range: 0.1...1.0, defaultValue: 0.4,
                               scale: .linear, decimalPlaces: 2)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "A：子弹（P₁₂ = P₁ + P₂，无相位结构）",
                xAxis: .init(label: "屏位置 x（自然单位）"),
                yAxis: .init(label: "概率密度"),
                seriesNames: ["子弹"])),
            .lineSeries(LineSeriesSpec(
                title: "B：两独立波源（I₁₂ = I₁ + I₂，时间平均）",
                xAxis: .init(label: "屏位置 x（自然单位）"),
                yAxis: .init(label: "概率密度"),
                seriesNames: ["水波"])),
            .lineSeries(LineSeriesSpec(
                title: "C：电子（P₁₂ = |φ₁ + φ₂|²，干涉）",
                xAxis: .init(label: "屏位置 x（自然单位）"),
                yAxis: .init(label: "概率密度"),
                seriesNames: ["电子"])),
            .lineSeries(LineSeriesSpec(
                title: "D：被监视电子（P₁₂ = P₁ + P₂，路径信息）",
                xAxis: .init(label: "屏位置 x（自然单位）"),
                yAxis: .init(label: "概率密度"),
                seriesNames: ["被监视电子"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let d = input.slider("d")
        let a = input.slider("a")
        let (x, bullets, waveB, electron, watched) =
            QuantumBehaviorMath.distributions(d: d, a: a)

        let fringe = 1.0 / d              // 条纹间距 λL/d（λ=L=1）
        let firstMin = 0.5 / d            // 一级极小 λL/(2d)

        func chart(_ name: String, _ ys: [Double], refs: [ReferenceLine] = []) -> ChartData {
            .lineSeries(LineSeriesData(
                spec: LineSeriesSpec(
                    xAxis: .init(label: "屏位置 x（自然单位）"),
                    yAxis: .init(label: "概率密度"),
                    seriesNames: [name]),
                series: [.init(name: name, points: Num.strided(x, ys, stride: 6))],
                referenceLines: refs))
        }

        // 50 kV 电子德布罗意波长（脚本附加输出：h / √(2·m_e·e·V)）
        let h = try constants.value("h")
        let me = try constants.value("m_e")
        let e = try constants.value("e")
        let lam50kV = h / (2 * me * e * 50e3).squareRoot()

        return SimResult(
            charts: [
                chart("子弹", bullets),
                chart("水波", waveB),
                chart("电子", electron, refs: [ReferenceLine(
                    id: "fmin", label: String(format: "一级极小 x = +%.3f", firstMin),
                    axis: .x, value: firstMin)]),
                chart("被监视电子", watched),
            ],
            summary: [
                .init(id: "fringe", title: "条纹间距 Δx",
                      value: String(format: "%.3f λL/d", fringe), note: "λL/d"),
                .init(id: "fmin", title: "一级极小",
                      value: String(format: "x = ±%.3f", firstMin), note: "λL/(2d)"),
                .init(id: "deb", title: "λ(电子, 50 kV)",
                      value: String(format: "%.2f pm", lam50kV * 1e12),
                      note: "h/√(2m_e·eV)"),
            ],
            theory: TheoryCard(
                title: "量子行为三实验",
                formulas: [
                    "子弹：P₁₂ = P₁ + P₂（概率直接相加，无相位）",
                    "电子：P₁₂ = |φ₁ + φ₂|² = 4·sinc²(πax/λL)·cos²(πdx/λL)",
                    "监视：路径信息获取 → 相位破坏 → 回到 P₁ + P₂",
                    "德布罗意波长 λ = h/p = h/√(2m·eV)",
                ],
                reading: "对比 C 与 D：同样的几何，是否获取路径信息决定干涉条纹存灭——观测本身改变结果。"))
    }
}
