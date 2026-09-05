import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/腔量子电动力学_JaynesCummings.py 与 _量子比特Rabi.py）

/// 腔 QED 纯函数核：共振真空 Rabi 振荡、缀饰态吸收谱、受驱量子比特退相干振荡。
enum CavityQEDMath {

    /// 共振真空 Rabi 振荡激发概率 P_e(t) = cos²(g t) = (1 + cos Ω_R t)/2（n = 0 光子）。
    static func vacuumRabiPe(_ t: Double, g: Double) -> Double {
        let c = cos(g * t)
        return c * c
    }

    /// n 光子 JC Rabi 频率 Ω_n = 2g√(n+1)（真空 n=0 → Ω_R = 2g）。
    static func rabiFrequency(g: Double, n: Int = 0) -> Double {
        2.0 * g * Double(n + 1).squareRoot()
    }

    /// 真空 Rabi 振荡周期 T_R = 2π/Ω_R = π/g。
    static func oscillationPeriod(g: Double) -> Double {
        .pi / g
    }

    /// 缀饰态吸收谱（双洛伦兹）：A(δ) = g²/((δ−g)²+(κ/2)²) + g²/((δ+g)²+(κ/2)²)。
    static func absorption(_ delta: Double, g: Double, kappa: Double) -> Double {
        let half = kappa / 2.0
        let d1 = delta - g, d2 = delta + g
        return g * g / (d1 * d1 + half * half) + g * g / (d2 * d2 + half * half)
    }

    /// 受驱量子比特 Rabi 振荡 P_e(t) = e^(−t/T₂)·sin²(Ω_R t/2)。
    static func drivenQubitPe(_ t: Double, omegaR: Double, t2: Double) -> Double {
        let s = sin(omegaR * t / 2.0)
        return exp(-t / t2) * s * s
    }

    /// 相干时间内可完成的振荡次数 N_osc ≈ Ω_R·T₂/π。
    static func coherentOscillationCount(omegaR: Double, t2: Double) -> Double {
        omegaR * t2 / .pi
    }

    /// 腔 QED（里德伯原子）与电路 QED（transmon）代表耦合（脚本固定值）。
    static let fgCQED = 47e3          // Hz，g/2π
    static let fgCircuitDefault = 5e6 // Hz，g/2π
}

// MARK: - 模块一：Jaynes-Cummings

/// 笔记 36《腔量子电动力学与量子比特实现》：共振真空 Rabi 振荡
/// P_e(t) = cos²(gt) 与缀饰态吸收谱双峰分裂（峰距 2g = Ω_R）。
struct JaynesCummingsModule: SimModule {

    let meta = ModuleMeta(
        id: "腔量子电动力学_JaynesCummings", title: "Jaynes-Cummings · 真空 Rabi 振荡",
        subtitle: "P_e(t)=cos²(gt) 与缀饰态吸收谱双峰分裂 2g——强耦合标志",
        category: .quantumOptics, noteNumber: 36, tier: .realtime, difficulty: .basic,
        keywords: ["Jaynes-Cummings", "真空 Rabi 振荡", "缀饰态", "强耦合", "腔量子电动力学",
                   "吸收谱", "vacuum Rabi", "dressed states", "normal mode splitting"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "fg", title: "电路 QED 耦合", symbol: "g/2π", unit: "MHz",
                               range: 0.1...50, defaultValue: 5,
                               scale: .log, decimalPlaces: 2)),
            .slider(SliderSpec(key: "kappa", title: "腔衰减率", symbol: "κ", unit: "MHz",
                               range: 0.01...5, defaultValue: 0.5,
                               scale: .log, decimalPlaces: 2)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "共振真空 Rabi 振荡 P_e(t)=cos²(gt)",
                xAxis: .init(label: "时间 t (μs)"),
                yAxis: .init(label: "激发概率 P_e(t)"),
                seriesNames: ["腔 QED（里德伯，g/2π=47 kHz）", "电路 QED（transmon）"])),
            .lineSeries(LineSeriesSpec(
                title: "缀饰态吸收谱（归一化，峰距 = 2g）",
                xAxis: .init(label: "探针失谐 δ (MHz)"),
                yAxis: .init(label: "归一化吸收 A/A_max"),
                seriesNames: ["电路 QED", "腔 QED（kHz 失谐，×10⁻³ MHz）"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let fgCircuit = input.slider("fg") * 1e6        // Hz
        let kappaMHz = input.slider("kappa")

        let gCQED = 2.0 * .pi * CavityQEDMath.fgCQED
        let gCircuit = 2.0 * .pi * fgCircuit
        let omegaRCQED = 2.0 * gCQED
        let omegaRCircuit = 2.0 * gCircuit

        // 图 1a：各平台画 3 个完整周期（脚本 linspace(0, 6π/Ω_R, 600) 抽稀至 512）
        let tCQED = Num.linspace(0, 6.0 * .pi / omegaRCQED, count: 600)
        let tCircuit = Num.linspace(0, 6.0 * .pi / omegaRCircuit, count: 600)
        let peCQED = tCQED.map { CavityQEDMath.vacuumRabiPe($0, g: gCQED) }
        let peCircuit = tCircuit.map { CavityQEDMath.vacuumRabiPe($0, g: gCircuit) }

        // 图 1b：吸收谱（脚本 κ_circuit = 2π·0.5 MHz，腔 QED 用 2π·5 kHz）
        let kappaCircuit = 2.0 * .pi * kappaMHz * 1e6
        let kappaCQED = 2.0 * .pi * 5e3
        let deltaCircuit = Num.linspace(-1.6 * gCircuit, 1.6 * gCircuit, count: 800)
        let deltaCQED = Num.linspace(-1.6 * gCQED, 1.6 * gCQED, count: 800)
        var aCircuit = deltaCircuit.map { CavityQEDMath.absorption($0, g: gCircuit, kappa: kappaCircuit) }
        var aCQED = deltaCQED.map { CavityQEDMath.absorption($0, g: gCQED, kappa: kappaCQED) }
        let maxCircuit = aCircuit.max() ?? 1
        let maxCQED = aCQED.max() ?? 1
        aCircuit = aCircuit.map { $0 / maxCircuit }
        aCQED = aCQED.map { $0 / maxCQED }

        let tRCQEDus = CavityQEDMath.oscillationPeriod(g: gCQED) * 1e6
        let tRCircuitus = CavityQEDMath.oscillationPeriod(g: gCircuit) * 1e6
        let splitMHz = 2.0 * fgCircuit / 1e6

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: charts[0].lineSeriesSpec!,
                    series: [
                        .init(name: "腔 QED（里德伯，g/2π=47 kHz）",
                              points: Num.strided(tCQED.map { $0 * 1e6 }, peCQED, stride: 1)),
                        .init(name: "电路 QED（transmon）",
                              points: Num.strided(tCircuit.map { $0 * 1e6 }, peCircuit, stride: 1),
                              colorIndex: 1),
                    ],
                    referenceLines: [
                        .init(label: "P_e = 1（初始激发）", axis: .y, value: 1, style: .subtle),
                    ])),
                .lineSeries(LineSeriesData(
                    spec: charts[1].lineSeriesSpec!,
                    series: [
                        .init(name: "电路 QED",
                              points: Num.strided(deltaCircuit.map { $0 / (2 * .pi * 1e6) }, aCircuit, stride: 1),
                              colorIndex: 1),
                        .init(name: "腔 QED（kHz 失谐，×10⁻³ MHz）",
                              points: Num.strided(deltaCQED.map { $0 / (2 * .pi * 1e3) * 1e-3 }, aCQED, stride: 1)),
                    ],
                    referenceLines: [
                        .init(label: "共振 δ = 0", axis: .x, value: 0, style: .subtle),
                        .init(label: "+g = \(String(format: "%.2f", fgCircuit / 1e6)) MHz",
                              axis: .x, value: fgCircuit / 1e6, style: .subtle),
                        .init(label: "−g = \(String(format: "%.2f", fgCircuit / 1e6)) MHz",
                              axis: .x, value: -fgCircuit / 1e6, style: .subtle),
                    ])),
            ],
            summary: [
                .init(id: "split", title: "缀饰态分裂 2g/2π",
                      value: String(format: "%.2f MHz", splitMHz),
                      note: "吸收谱双峰间距 = 真空 Rabi 频率"),
                .init(id: "TRcircuit", title: "电路 QED 振荡周期 T_R",
                      value: String(format: "%.3f μs", tRCircuitus),
                      note: "π/g（g/2π = \(String(format: "%.2f", fgCircuit / 1e6)) MHz）"),
                .init(id: "TRcqed", title: "腔 QED 振荡周期 T_R",
                      value: String(format: "%.3f μs", tRCQEDus),
                      note: "π/g（g/2π = 47 kHz，里德伯原子代表值）"),
                .init(id: "regime", title: "耦合区判定",
                      value: kappaMHz < fgCircuit ? "强耦合 g > κ" : "弱耦合 g ≤ κ",
                      note: "g/κ = \(String(format: "%.2f", fgCircuit / kappaMHz))（强耦合下方可分辨双峰）"),
            ],
            theory: TheoryCard(
                title: "Jaynes-Cummings 模型（笔记 36 式 1-5、7）",
                formulas: [
                    "共振真空 Rabi 振荡：P_e(t) = cos²(g t) = (1 + cos Ω_R t)/2，Ω_R = 2g",
                    "n 光子推广：Ω_n = 2g√(n+1)（真空 n=0 即上式）",
                    "缀饰态吸收谱：A(δ) = g²/((δ−g)²+(κ/2)²) + g²/((δ+g)²+(κ/2)²)",
                    "双峰间距 = 2g = Ω_R：强耦合区（g ≫ κ）的正常模分裂",
                ],
                reading: "左图：真空场中原子-腔交换能量，P_e 以 Ω_R = 2g 振荡——电路 QED（红）"
                    + "比腔 QED（蓝）快 100 倍。右图：弱探针扫描给出双洛伦兹吸收峰，"
                    + "峰距恰为 2g——单光子耦合可直接从谱上读出。"))
    }
}

// MARK: - 模块二：受驱量子比特 Rabi

/// 笔记 36《腔量子电动力学与量子比特实现》：受驱量子比特 Rabi 振荡
/// P_e(t) = e^(−t/T₂)·sin²(Ω_R t/2)——三平台相干时间尺度对比。
struct QubitRabiModule: SimModule {

    /// 平台相干时间（脚本代表值）。
    struct Platform {
        let name: String
        let t2: Double       // s
        let tMax: Double     // s（画图窗口）
        let colorIndex: Int
    }

    static let platforms: [Platform] = [
        .init(name: "离子阱（¹⁷¹Yb⁺/⁴⁰Ca⁺，T₂=1 s）", t2: 1.0, tMax: 4e-3, colorIndex: 0),
        .init(name: "超导 transmon（T₂=50 μs）", t2: 50e-6, tMax: 200e-6, colorIndex: 1),
        .init(name: "NV 色心（室温，T₂=1 ms）", t2: 1e-3, tMax: 4e-3, colorIndex: 2),
    ]

    let meta = ModuleMeta(
        id: "腔量子电动力学_量子比特Rabi", title: "量子比特 Rabi · 退相干包络",
        subtitle: "P_e(t)=e^(−t/T₂)·sin²(Ω_R t/2)——三平台相干时间尺度对比",
        category: .quantumOptics, noteNumber: 36, tier: .realtime, difficulty: .basic,
        keywords: ["量子比特", "Rabi 振荡", "退相干", "相干时间", "T₂", "离子阱", "transmon", "NV 色心",
                   "qubit", "decoherence", "coherence time"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "fdrive", title: "驱动频率", symbol: "Ω_R/2π", unit: "MHz",
                               range: 0.1...50, defaultValue: 1,
                               scale: .log, decimalPlaces: 2)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "受驱量子比特 Rabi 振荡（含退相干包络）",
                xAxis: .init(label: "时间 t (μs)"),
                yAxis: .init(label: "P_e(t) = e^(−t/T₂)·sin²(Ω_R t/2)"),
                seriesNames: QubitRabiModule.platforms.map(\.name))),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let fDrive = input.slider("fdrive") * 1e6     // Hz
        let omegaR = 2.0 * .pi * fDrive

        var series: [SeriesPoints] = []
        for p in QubitRabiModule.platforms {
            let t = Num.linspace(0, p.tMax, count: 2000)
            let pe = t.map { CavityQEDMath.drivenQubitPe($0, omegaR: omegaR, t2: p.t2) }
            series.append(.init(name: p.name,
                                points: Num.strided(t.map { $0 * 1e6 }, pe, stride: 4),
                                colorIndex: p.colorIndex))
        }

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: charts[0].lineSeriesSpec!,
                    series: series,
                    referenceLines: [
                        .init(label: "P_e = 1（无退相干上限）", axis: .y, value: 1, style: .subtle),
                    ])),
            ],
            summary: [
                .init(id: "osc", title: "相干时间内振荡次数",
                      value: String(format: "%.2e", CavityQEDMath.coherentOscillationCount(
                          omegaR: omegaR, t2: QubitRabiModule.platforms[0].t2)),
                      note: "离子阱 T₂=1 s：Ω_R·T₂/π"),
                .init(id: "period", title: "Rabi 振荡周期",
                      value: String(format: "%.3f μs", 2.0 * .pi / omegaR * 1e6),
                      note: "2π/Ω_R（Ω_R/2π = \(String(format: "%.2f", fDrive / 1e6)) MHz）"),
                .init(id: "T2ion", title: "离子阱 T₂",
                      value: "1 s 量级", note: "¹⁷¹Yb⁺ / ⁴⁰Ca⁺（可至 1-10 s）"),
                .init(id: "T2sc", title: "超导 transmon T₂",
                      value: "50 μs 量级", note: "10-100 μs（探测优化后可延长）"),
                .init(id: "T2nv", title: "NV 色心 T₂",
                      value: "1 ms 量级", note: "室温；低温同位素纯化后可达 ~s"),
            ],
            theory: TheoryCard(
                title: "受驱量子比特 Rabi 振荡（笔记 36 式 10-11）",
                formulas: [
                    "P_e(t) = e^(−t/T₂)·sin²(Ω_R t/2)——相干振荡 × 退相干包络",
                    "包络 e^(−t/T₂)：T₂ 为横向相干时间",
                    "相干振荡次数 N_osc ≈ Ω_R·T₂/π",
                    "平台量级：离子阱 ~1 s；transmon ~50 μs；NV 室温 ~1 ms",
                ],
                reading: "三条曲线共享同一驱动频率，差别只在 T₂：离子阱（绿）窗口 4 ms 内"
                    + "维持上万次振荡；transmon（红）窗口 200 μs 内几十次振荡即被包络压平；"
                    + "NV（紫）居中。量子门必须在 T₂ 内完成——相干时间即「计算预算」。"))
    }
}
