import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/角动量与自旋_斯特恩-盖拉赫分裂.py）

/// 斯特恩-盖拉赫空间量子化偏转的纯函数核。
enum SternGerlachMath {

    /// 银原子质量（107.8682 u，Ag-107/109 天然平均）。
    static func silverMass(u: Double) -> Double { 107.8682 * u }

    /// 炉温热速度 v = √(2kBT/m)。
    static func thermalSpeed(T: Double, kB: Double, mass: Double) -> Double {
        (2 * kB * T / mass).squareRoot()
    }

    /// 偏转链：a = μB·∇B/m → 磁内 d₁ = ½at₁² → 漂移 d₂ = v_z·t₂ → dz = d₁+d₂。
    static func deflection(
        v: Double, muB: Double, gradB: Double, mass: Double,
        magnetLength L: Double, drift D: Double
    ) -> (a: Double, d1: Double, d2: Double, dz: Double) {
        let a = muB * gradB / mass
        let t1 = L / v
        let vz = a * t1
        let d1 = 0.5 * a * t1 * t1
        let d2 = vz * (D / v)
        return (a, d1, d2, d1 + d2)
    }

    /// 单束高斯：I±(z) = exp(−((z∓dz)/(√2σ))²)，σ = 0.25 dz。
    static func beamIntensity(z: Double, dz: Double, up: Bool) -> Double {
        let sigma = 0.25 * dz
        let d = z - (up ? dz : -dz)
        let t = d / (2.0.squareRoot() * sigma)
        return exp(-t * t)
    }

    /// 探测屏网格：z = linspace(−3dz, 3dz, 600)。
    static func detectorGrid(dz: Double) -> [Double] {
        Num.linspace(-3 * dz, 3 * dz, count: 600)
    }
}

// MARK: - 模块

/// 笔记 16《角动量与自旋》§6.2：银原子束在不均匀磁场中的空间量子化——
/// m_s = ±1/2 两束分立偏转（非经典连续），量级与 1922 年实验自洽。
struct SternGerlachModule: SimModule {

    let meta = ModuleMeta(
        id: "角动量与自旋_斯特恩-盖拉赫分裂", title: "斯特恩-盖拉赫 · 空间量子化",
        subtitle: "自旋 1/2 的银原子束被梯度磁场劈成两束——角动量分立的直接证据",
        category: .angularCentral, noteNumber: 16, tier: .seconds, difficulty: .basic,
        keywords: ["斯特恩-盖拉赫", "空间量子化", "自旋", "磁矩", "Stern-Gerlach",
                   "m_s", "梯度磁场", "银原子"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "T", title: "炉温", symbol: "T", unit: "K",
                               range: 500...1500, defaultValue: 1000,
                               decimalPlaces: 0)),
            .slider(SliderSpec(key: "gradB", title: "磁场梯度", symbol: "dB/dz",
                               unit: "T/m", range: 100...10000, defaultValue: 1000,
                               scale: .log, decimalPlaces: 0)),
            .slider(SliderSpec(key: "Lm", title: "磁铁长度", symbol: "L", unit: "m",
                               range: 0.05...0.30, defaultValue: 0.10,
                               decimalPlaces: 2)),
            .slider(SliderSpec(key: "D", title: "漂移距离", symbol: "D", unit: "m",
                               range: 0.2...2.0, defaultValue: 1.0,
                               decimalPlaces: 1)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "探测屏强度分布（两束分立）",
                xAxis: .init(label: "探测器位置 z (mm)"), yAxis: .init(label: "相对束强"),
                seriesNames: ["m_s = +1/2", "m_s = −1/2", "总强度"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let T = input.slider("T")
        let gradB = input.slider("gradB")
        let Lm = input.slider("Lm")
        let D = input.slider("D")
        let u = try constants.value("u"), kB = try constants.value("kB")
        let muB = try constants.value("mu_B")

        let mAg = SternGerlachMath.silverMass(u: u)
        let v = SternGerlachMath.thermalSpeed(T: T, kB: kB, mass: mAg)
        let (a, d1, d2, dz) = SternGerlachMath.deflection(
            v: v, muB: muB, gradB: gradB, mass: mAg, magnetLength: Lm, drift: D)

        // 探测屏曲线（脚本原样 600 网格；显示抽稀 stride 3 → 200 点）
        let grid = SternGerlachMath.detectorGrid(dz: dz)
        var upPts = [Point](), dnPts = [Point](), totPts = [Point]()
        upPts.reserveCapacity(200); dnPts.reserveCapacity(200); totPts.reserveCapacity(200)
        for i in grid.indices where i % 3 == 0 {
            let iu = SternGerlachMath.beamIntensity(z: grid[i], dz: dz, up: true)
            let idn = SternGerlachMath.beamIntensity(z: grid[i], dz: dz, up: false)
            upPts.append(Point(x: grid[i] * 1e3, y: iu))
            dnPts.append(Point(x: grid[i] * 1e3, y: idn))
            totPts.append(Point(x: grid[i] * 1e3, y: iu + idn))
        }

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: LineSeriesSpec(
                        title: "探测屏强度分布（两束分立）",
                        xAxis: .init(label: "探测器位置 z (mm)"), yAxis: .init(label: "相对束强"),
                        seriesNames: ["m_s = +1/2", "m_s = −1/2", "总强度"]),
                    series: [
                        .init(name: "m_s = +1/2", points: upPts),
                        .init(name: "m_s = −1/2", points: dnPts, colorIndex: 1),
                        .init(name: "总强度", points: totPts, colorIndex: 2),
                    ])),
            ],
            summary: [
                .init(id: "v", title: "热速度 v",
                      value: String(format: "%.1f", v), note: "m/s（炉温决定）"),
                .init(id: "a", title: "横向加速度 a",
                      value: String(format: "%.3e", a), note: "m/s²（F = μz·∇B）"),
                .init(id: "dz", title: "单束偏移 dz",
                      value: String(format: "%.3f", dz * 1e3), note: "mm（磁内 "
                          + String(format: "%.3f", d1 * 1e3) + " + 漂移 "
                          + String(format: "%.3f", d2 * 1e3) + "）"),
                .init(id: "sep", title: "两束分离 2dz",
                      value: String(format: "%.3f", 2 * dz * 1e3), note: "mm"),
            ],
            theory: TheoryCard(
                title: "斯特恩-盖拉赫：自旋的发现",
                formulas: [
                    "F = μz·dB/dz：磁矩在梯度场受力（Ag 价电子 L = 0，μ = −gs μB S）",
                    "dz = ½·a·(L/v)² + a·(L/v)·(D/v)：磁内抛物 + 漂移直线",
                    "m_s = ±1/2 ⇒ 屏上两束（经典轨道角动量应为连续环）",
                    "两束等强 ⇒ 自旋态 |↑⟩/|↓⟩ 各占 1/2（未极化束）",
                ],
                reading: "1922 年的银原子实验：屏上两条分立沉积，宣告角动量空间量子化。"
                    + "拖大 dB/dz 或 D 看分离拉大；升温 T 让束更快、偏转反而变小——"
                    + "偏转 ∝ 1/v² 量级，热速度是分离的敌人。"))
    }
}
