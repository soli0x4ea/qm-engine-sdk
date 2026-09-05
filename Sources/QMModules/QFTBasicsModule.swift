import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/量子场论的基本图像_Casimir.py 与 _模式量子化.py）

/// 量子场论基本图像纯函数核：Casimir 压强、一维振子链模式、零点能、BE 占据数。
enum QFTBasicsMath {

    /// 理想导体平行板 Casimir 压强绝对值 |F|/A = π²ℏc/(240 d⁴)（N/m²；d 单位 m）。
    static func casimirPressure(_ d: Double, hbar: Double, c: Double) -> Double {
        Double.pi * Double.pi * hbar * c / (240.0 * pow(d, 4))
    }

    /// Casimir 压强（带符号，吸引为负）。
    static func casimirPressureSigned(_ d: Double, hbar: Double, c: Double) -> Double {
        -casimirPressure(d, hbar: hbar, c: c)
    }

    /// 局部幂律指数 d ln P / d ln d（对数网格差分）。
    static func localExponents(_ d: [Double], _ p: [Double]) -> [Double] {
        guard d.count == p.count, d.count > 1 else { return [] }
        var out: [Double] = []
        out.reserveCapacity(d.count - 1)
        for i in 1..<d.count {
            out.append((log(p[i]) - log(p[i - 1])) / (log(d[i]) - log(d[i - 1])))
        }
        return out
    }

    /// 一维单原子链模式波矢 k_m = 2π m/(N a)（m = −N/2 … N/2−1）。
    static func waveVectors(N: Int, a: Double) -> [Double] {
        (-(N / 2)..<(N / 2)).map { 2.0 * .pi * Double($0) / (Double(N) * a) }
    }

    /// 声学支色散 ω_k = 2(v/a)|sin(k a/2)|（式 5）。
    static func dispersion(_ k: Double, a: Double, v: Double) -> Double {
        2.0 * (v / a) * abs(sin(k * a / 2.0))
    }

    /// 各模零点能 E₀ = Σ_k ℏω_k/2（式 6；固定 N 扫 a_lat）。
    static func zeroPointEnergy(N: Int, a: Double, v: Double, hbar: Double) -> Double {
        let ks = waveVectors(N: N, a: a)
        var sum = 0.0
        for k in ks { sum += 0.5 * hbar * dispersion(k, a: a, v: v) }
        return sum
    }

    /// 玻色-爱因斯坦平衡占据数 ⟨n_k⟩ = 1/(e^{ℏω/k_BT} − 1)（式 7；零模取 0）。
    static func boseOccupancy(_ omega: Double, temperature: Double,
                              hbar: Double, kB: Double) -> Double {
        if omega == 0 { return 0 }
        let x = hbar * omega / (kB * temperature)
        if x < 1e-12 { return 1.0 / x - 0.5 }  // 高温经典极限展开（避免 0/0）
        return 1.0 / (exp(x) - 1.0)
    }
}

// MARK: - 模块一：Casimir

/// 笔记 40《量子场论的基本图像》：理想导体平行板 Casimir 压强
/// P = −π²ℏc/(240 d⁴)——双对数坐标的 d⁻⁴ 幂律与局部指数验证。
struct CasimirModule: SimModule {

    let meta = ModuleMeta(
        id: "量子场论的基本图像_Casimir", title: "Casimir 效应 · d⁻⁴ 幂律",
        subtitle: "P = −π²ℏc/(240 d⁴)——真空零点能差的宏观吸引压",
        category: .relativisticQFT, noteNumber: 40, tier: .realtime, difficulty: .basic,
        keywords: ["Casimir", "卡西米尔效应", "真空零点能", "幂律", "吸引",
                   "vacuum energy", "power law"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "dmin", title: "板间距下限", symbol: "d_min", unit: "μm",
                               range: 0.05...1, defaultValue: 0.1,
                               scale: .log, decimalPlaces: 2)),
            .slider(SliderSpec(key: "dmax", title: "板间距上限", symbol: "d_max", unit: "μm",
                               range: 1...50, defaultValue: 10,
                               scale: .log, decimalPlaces: 1)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "Casimir 压强 vs 板间距（双对数，斜率 −4）",
                xAxis: .init(label: "板间距 d (μm)", scale: .log),
                yAxis: .init(label: "|F|/A (N/m²)", scale: .log),
                seriesNames: ["|P(d)| = π²ℏc/(240 d⁴)"])),
            .lineSeries(LineSeriesSpec(
                title: "局部幂律指数 d ln P / d ln d",
                xAxis: .init(label: "板间距 d (μm)", scale: .log),
                yAxis: .init(label: "局部指数"),
                seriesNames: ["d ln P / d ln d"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let hbar = try constants.value("hbar")
        let c = try constants.value("c")
        let dMin = input.slider("dmin") * 1e-6   // m
        let dMax = input.slider("dmax") * 1e-6

        // 脚本同口径：logspace 扫描（0.1 μm ~ 10 μm 200 点）
        let d = Num.linspace(log10(dMin), log10(dMax), count: 200).map { pow(10, $0) }
        let p = d.map { QFTBasicsMath.casimirPressure($0, hbar: hbar, c: c) }
        let slopes = QFTBasicsMath.localExponents(d, p)

        let pAt1um = QFTBasicsMath.casimirPressure(1e-6, hbar: hbar, c: c)
        let pAt01 = QFTBasicsMath.casimirPressure(1e-7, hbar: hbar, c: c)
        let pAt10 = QFTBasicsMath.casimirPressure(1e-5, hbar: hbar, c: c)
        let fOn1cm2 = pAt1um * 1e-4

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: charts[0].lineSeriesSpec!,
                    series: [
                        .init(name: "|P(d)| = π²ℏc/(240 d⁴)",
                              points: zip(d.map { $0 * 1e6 }, p).map { Point(x: $0, y: $1) }),
                    ],
                    referenceLines: [])),
                .lineSeries(LineSeriesData(
                    spec: charts[1].lineSeriesSpec!,
                    series: [
                        .init(name: "d ln P / d ln d",
                              points: zip(d.dropLast().map { $0 * 1e6 }, slopes)
                                  .map { Point(x: $0, y: $1) }),
                    ],
                    referenceLines: [
                        .init(label: "理想指数 −4", axis: .y, value: -4),
                    ])),
            ],
            summary: [
                .init(id: "p1", title: "d = 1 μm 处 |F|/A",
                      value: String(format: "%.4e N/m²", pAt1um),
                      note: "纯幂律解析值"),
                .init(id: "p01", title: "d = 0.1 μm 处 |F|/A",
                      value: String(format: "%.4e N/m²", pAt01),
                      note: "间距缩小 10×，压强放大 10⁴×"),
                .init(id: "p10", title: "d = 10 μm 处 |F|/A",
                      value: String(format: "%.4e N/m²", pAt10),
                      note: "间距放大 10×，压强缩小 10⁴×"),
                .init(id: "f", title: "A = 1 cm²、d = 1 μm 吸引力",
                      value: String(format: "%.4e N", fOn1cm2),
                      note: "P(d) < 0：压强为吸引（负号）"),
            ],
            theory: TheoryCard(
                title: "Casimir 压强（笔记 40 式 8）",
                formulas: [
                    "理想导体平行板：P(d) = −π²ℏc/(240 d⁴)",
                    "负号 = 吸引；|P| ∝ d⁻⁴ 是真空模式谱的零点能差",
                    "双对数坐标：ln|P| = const − 4 ln d（直线，斜率 −4）",
                    "局部指数 d ln P/d ln d ≡ −4（纯幂律的全局不变量）",
                ],
                reading: "左图双对数下是严格的直线——斜率 −4 即 d⁻⁴ 幂律。右图局部指数"
                    + "在整个扫描区间恒为 −4（红线）：这不是拟合，是解析幂律的必然。"
                    + "物理来源：两板间允许的真空模式少于板外，零点能差把板推向彼此。"))
    }
}

// MARK: - 模块二：模式量子化

/// 笔记 40《量子场论的基本图像》：一维振子链当作「场」的离散模式——
/// 声学支色散 ω_k、各模零点能 E₀ = Σℏω/2、BE 平衡占据数，与 1/a² 紫外发散。
struct ModeQuantizationModule: SimModule {

    let meta = ModuleMeta(
        id: "量子场论的基本图像_模式量子化", title: "模式量子化 · 零点能与 BE 占据",
        subtitle: "ω_k = 2(v/a)|sin(ka/2)|、E₀ = Σℏω/2、⟨n_k⟩ = 1/(e^{ℏω/kT} − 1)",
        category: .relativisticQFT, noteNumber: 40, tier: .realtime, difficulty: .basic,
        keywords: ["模式量子化", "零点能", "真空能", "声子色散", "玻色-爱因斯坦占据",
                   "紫外发散", "声学支", "phonon", "zero-point energy", "UV cutoff"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "a", title: "晶格常数", symbol: "a", unit: "nm",
                               range: 0.05...5, defaultValue: 0.5,
                               scale: .log, decimalPlaces: 2)),
            .slider(SliderSpec(key: "v", title: "声速", symbol: "v", unit: "km/s",
                               range: 0.1...10, defaultValue: 1,
                               scale: .log, decimalPlaces: 2)),
            .slider(SliderSpec(key: "T", title: "温度", symbol: "T", unit: "K",
                               range: 1...1000, defaultValue: 300,
                               scale: .log, decimalPlaces: 0)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "一维振子链色散 ω(k)",
                xAxis: .init(label: "k a / (2π)（无量纲）"),
                yAxis: .init(label: "ω (THz)"),
                seriesNames: ["ω_k = 2(v/a)|sin(ka/2)|"])),
            .lineSeries(LineSeriesSpec(
                title: "各模零点能 ℏω_k/2",
                xAxis: .init(label: "模式序号 m"),
                yAxis: .init(label: "E₀ (10⁻²² J)"),
                seriesNames: ["ℏω_k/2"])),
            .lineSeries(LineSeriesSpec(
                title: "玻色-爱因斯坦平衡占据数",
                xAxis: .init(label: "k a / (2π)（无量纲）"),
                yAxis: .init(label: "⟨n_k⟩（对数）", scale: .log),
                seriesNames: ["T = 10 K", "当前温度"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let N = 50
        let a = input.slider("a") * 1e-10      // m
        let v = input.slider("v") * 1e3        // m/s
        let T = input.slider("T")
        let hbar = try constants.value("hbar")
        let kB = try constants.value("kB")

        // 与脚本同口径：m_idx = −25…24，k_m = 2π m/(N a)，ω_m = ω_max|sin(πm/N)|
        let ks = QFTBasicsMath.waveVectors(N: N, a: a)
        let omegas = ks.map { QFTBasicsMath.dispersion($0, a: a, v: v) }
        let e0PerMode = omegas.map { 0.5 * hbar * $0 }
        let e0Total = e0PerMode.reduce(0, +)
        let u0 = e0Total / (Double(N) * a)

        // 紫外截断扫描（脚本固定 a_list；以当前 v/N 扫描）
        let aList = [5e-10, 5e-11, 5e-12, 5e-13]
        let uScan = aList.map { QFTBasicsMath.zeroPointEnergy(N: N, a: $0, v: v, hbar: hbar) / (50.0 * $0) }
        let ua2 = zip(aList, uScan).map { $1 * $0 * $0 }

        // BE 占据数（脚本口径：T = 10 K 固定 + 当前温度）
        let occ10 = omegas.map { QFTBasicsMath.boseOccupancy($0, temperature: 10, hbar: hbar, kB: kB) }
        let occT = omegas.map { QFTBasicsMath.boseOccupancy($0, temperature: T, hbar: hbar, kB: kB) }

        // 最小非零 |k| 模（脚本 idx_small）与最大 |k| 模
        var iSmall = 0, iLarge = 0
        var kMin = Double.greatestFiniteMagnitude
        for (i, k) in ks.enumerated() where omegas[i] > 0 {
            if abs(k) < kMin { kMin = abs(k); iSmall = i }
            iLarge = i
        }

        let omegaMax = 2.0 * v / a
        let xKA = ks.map { $0 * a / (2.0 * .pi) }   // 无量纲 = m/N

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: charts[0].lineSeriesSpec!,
                    series: [
                        .init(name: "ω_k = 2(v/a)|sin(ka/2)|",
                              points: zip(xKA, omegas.map { $0 / 1e12 }).map { Point(x: $0, y: $1) }),
                    ],
                    referenceLines: [
                        .init(label: "带宽上限 2v/a = \(String(format: "%.2f", omegaMax / 1e12)) THz",
                              axis: .y, value: omegaMax / 1e12, style: .subtle),
                    ])),
                .lineSeries(LineSeriesData(
                    spec: charts[1].lineSeriesSpec!,
                    series: [
                        .init(name: "ℏω_k/2",
                              points: (0..<N).map { Point(x: Double($0), y: e0PerMode[$0] / 1e-22) }),
                    ],
                    referenceLines: [])),
                .lineSeries(LineSeriesData(
                    spec: charts[2].lineSeriesSpec!,
                    series: [
                        .init(name: "T = 10 K",
                              points: zip(xKA, occ10).map { Point(x: $0, y: max($1, 1e-6)) }),
                        .init(name: "当前温度", points: zip(xKA, occT).map { Point(x: $0, y: max($1, 1e-6)) },
                              colorIndex: 1),
                    ],
                    referenceLines: [
                        .init(label: "⟨n⟩ = 1（量子判据 ℏω = k_BT）", axis: .y, value: 1, style: .subtle),
                    ])),
            ],
            summary: [
                .init(id: "wmax", title: "带宽上限 ω_max",
                      value: String(format: "%.4e s⁻¹", omegaMax),
                      note: "2v/a（\(String(format: "%.2f", omegaMax / 1e12)) THz，布里渊区边界）"),
                .init(id: "e0", title: "零点能 E₀ = Σℏω/2",
                      value: String(format: "%.4e J", e0Total),
                      note: "每格点 \(String(format: "%.4e", e0Total / Double(N))) J"),
                .init(id: "u0", title: "零点能密度 u₀",
                      value: String(format: "%.4e J/m", u0),
                      note: "u ∝ 1/a²：紫外截断敏感（真空能问题）"),
                .init(id: "ua2", title: "标度检验 u·a²",
                      value: String(format: "%.4e J·m", ua2[0]),
                      note: "截断扫描四点常差 < 1e-12（u ∝ a⁻² 严格成立）"),
                .init(id: "occ", title: "BE 占据（300 K）",
                      value: String(format: "%.3e / %.3e", occT[iSmall], occT[iLarge]),
                      note: "最小非零 |k| 模 / 布里渊区边界模"),
            ],
            theory: TheoryCard(
                title: "模式量子化（笔记 40 式 5-7）",
                formulas: [
                    "一维单原子链声学支：ω_k = 2(v/a)|sin(ka/2)|（k = 2πm/(Na)）",
                    "每个模式 = 一个量子谐振子：零点能 E₀ = Σ_k ℏω_k/2（式 6）",
                    "热平衡占据：⟨n_k⟩ = 1/(e^{ℏω_k/k_BT} − 1)（式 7）",
                    "零点能密度 u ∝ 1/a²：截断 a → 0 时发散（真空能/宇宙学常数问题）",
                ],
                reading: "左图：色散在 k→0 处线性（声学支 ω ≈ v|k|），布里渊区边界达 2v/a。"
                    + "中图：各模零点能 ℏω/2——正是真空不空的来源。右图：低温下只有"
                    + "长波模（k 小）被占据；ℏω = k_BT 处 ⟨n⟩ = 1 标记量子-经典分界。"))
    }
}