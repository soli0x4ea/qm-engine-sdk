import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/反粒子与狄拉克海_对产生阈值.py 与 _狄拉克海洞.py）

/// 反粒子与狄拉克海纯函数核：对产生阈值（含核反冲修正）与狄拉克海洞示意。
enum DiracSeaMath {

    /// 单轻子静能 m_e c²（MeV）。
    static func restEnergy(constants: ConstantsSet) throws -> Double {
        try constants.value("me_c2_MeV")
    }

    /// 无反冲对产生阈值 E_th = 2 m_e c²（MeV）。
    static func pairThreshold(constants: ConstantsSet) throws -> Double {
        try 2.0 * restEnergy(constants: constants)
    }

    /// 含核反冲的严格阈值（实验室系、核静止）：
    /// s = M_N² + 2M_N·E_γ，阈值末态 (M_N + 2m_e)² ⇒ E_γ,min = 2 m_e c²(1 + m_e/M_N)。
    static func pairThresholdRecoil(constants: ConstantsSet, nucleusMass: Double) throws -> Double {
        let me = try constants.value("m_e")
        let c = try constants.value("c")
        let e = try constants.value("e")
        let restJ = me * c * c
        let restMeV = restJ / (1e6 * e)
        return 2.0 * restMeV * (1.0 + me / nucleusMass)
    }

    /// 狄拉克海占据态抽样（N 个点，均匀分布于负能连续谱示意区）。
    /// 种子 = 20260901（Python 脚本原值；SplitMix64 序列不同、分布相同，
    /// 规范见 docs/RNG_SEED_POLICY.md）。
    static func seaPoints(count n: Int, seed: UInt64 = 20_260_901,
                          eBound: Double = 1.30) -> [(x: Double, y: Double)] {
        var rng = SplitMix64(seed: seed)
        var pts: [(Double, Double)] = []
        pts.reserveCapacity(n)
        for _ in 0..<n {
            let x = -0.9 + 1.8 * rng.nextDouble()
            let y = -(eBound - 0.05) + (eBound - 0.10) * rng.nextDouble()
            pts.append((x, y))
        }
        return pts
    }
}

// MARK: - 模块一：对产生阈值

/// 笔记 39《反粒子与狄拉克海》：e⁺e⁻ 对产生阈值 E_th = 2 m_e c²
/// 及核库仑场三体阈值的反冲修正 E_γ,min = 2mc²(1 + m_e/M_N)——能轴示意。
struct PairProductionModule: SimModule {

    let meta = ModuleMeta(
        id: "反粒子与狄拉克海_对产生阈值", title: "对产生阈值 · 2mc²",
        subtitle: "E_th = 2 m_e c² ≈ 1.022 MeV——核场中单光子产阈的严格三体阈值含反冲修正",
        category: .relativisticQFT, noteNumber: 39, tier: .realtime, difficulty: .basic,
        keywords: ["对产生", "阈值", "正负电子对", "反冲修正", "核库仑场",
                   "pair production", "threshold", "recoil"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "A", title: "靶核质量数", symbol: "A", unit: "",
                               range: 1...250, defaultValue: 1,
                               step: 1, decimalPlaces: 0)),
            .slider(SliderSpec(key: "Eg", title: "入射光子能量", symbol: "E_γ", unit: "MeV",
                               range: 0...1.6, defaultValue: 1.2,
                               decimalPlaces: 3)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .schematic(SchematicSpec(
                title: "e⁺e⁻ 对产生阈值示意图",
                xAxis: .init(label: "（示意横轴）"),
                yAxis: .init(label: "光子 / 电子对能量 E (MeV)")))
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let A = input.slider("A")
        let eGamma = input.slider("Eg")
        let mP = try constants.value("m_p")
        let eRest = try DiracSeaMath.restEnergy(constants: constants)
        let eTh = try DiracSeaMath.pairThreshold(constants: constants)
        let eThRecoil = try DiracSeaMath.pairThresholdRecoil(
            constants: constants, nucleusMass: A * mP)
        let recoil = eThRecoil - eTh

        // W13（B4）：动态元素——入射光子标记随 E_γ 滑杆上下移动，
        // 跨越阈值线即「亚阈值 / 可产对」判定可视化。
        let feasible = eGamma >= eThRecoil
        let verdict = feasible ? "E_γ ≥ 含反冲阈 → 对产生可行"
                               : "E_γ < 含反冲阈 → 亚阈值（真空中无自由对产生）"

        return SimResult(
            charts: [
                .schematic(SchematicData(
                    spec: charts[0].schematicSpec!,
                    bands: [
                        .init(label: "亚阈值区（真空中无自由对产生）",
                              xRange: -1.2...1.2, yRange: 0...eTh,
                              filled: true, colorIndex: 4),
                    ],
                    markers: [
                        .init(label: "m_e c²（单轻子静能）", x: 0, y: eRest,
                              filled: true, colorIndex: 2),
                        // 入射光子（随 E_γ 滑杆移动的动态标记）
                        .init(label: String(format: "入射 γ：E_γ = %.3f MeV", eGamma),
                              x: -0.9, y: eGamma, filled: true, colorIndex: 3),
                    ],
                    callouts: [
                        .init(text: "2 m_e c² = \(String(format: "%.4f", eTh)) MeV（阈值）",
                              anchor: Point(x: 0.3, y: eTh),
                              arrowEnd: Point(x: 0.05, y: eTh), colorIndex: 0),
                        .init(text: "含反冲阈 \(String(format: "%.6f", eThRecoil)) MeV\n（A = \(Int(A))，修正 +\(String(format: "%.3e", recoil)) MeV）",
                              anchor: Point(x: 0.3, y: eThRecoil),
                              arrowEnd: Point(x: 0.05, y: eThRecoil), colorIndex: 1),
                        .init(text: verdict,
                              anchor: Point(x: 0, y: eGamma + 0.12),
                              arrowEnd: Point(x: -0.88, y: eGamma), colorIndex: 3),
                    ],
                    referenceLines: [
                        .init(label: "E = 0（静止）", axis: .y, value: 0, style: .subtle),
                    ])),
            ],
            summary: [
                .init(id: "eRest", title: "电子静能 m_e c²",
                      value: String(format: "%.8f MeV", eRest),
                      note: "CODATA \(constants.version.rawValue)"),
                .init(id: "eTh", title: "对产生阈值 2 m_e c²",
                      value: String(format: "%.8f MeV", eTh),
                      note: "无反冲近似（重靶 / 二体近似）"),
                .init(id: "recoil", title: "反冲修正 E_γ,min − 2mc²",
                      value: String(format: "%.6e MeV", recoil),
                      note: "2mc²·(m_e/M_N) > 0：重核靶修正趋零"),
                .init(id: "verdict", title: "当前 E_γ 判定",
                      value: feasible ? "可行（越阈）" : "不可行（亚阈值）",
                      note: String(format: "E_γ = %.3f MeV vs 含反冲阈 %.6f MeV", eGamma, eThRecoil)),
            ],
            theory: TheoryCard(
                title: "对产生阈值（笔记 39 第六节 a）",
                formulas: [
                    "二体阈值：E_th = 2 m_e c² ≈ 1.022 MeV",
                    "三体严格阈值（核静止实验室系）：s = M_N² + 2 M_N E_γ",
                    "阈值末态为单一静止系统 (M_N + 2 m_e)²",
                    "⇒ E_γ,min = 2 m_e c² (1 + m_e/M_N)——反冲修正恒正",
                ],
                reading: "红色带为亚阈值区：能量低于 2mc² 时真空中无法产生有质量对。"
                    + "阈值线之上单光子可在核库仑场中转化为 e⁺e⁻ 对（核吸收反冲动量）。"
                    + "换不同质量数 A 的靶：反冲修正 ∝ 1/A，重核极限回到干净的 2mc²。"))
    }
}

// MARK: - 模块二：狄拉克海洞

/// 笔记 39《反粒子与狄拉克海》：负能海填满 + 洞 = 正电子——
/// 电荷/能量符号对应表（洞 Q = +e、E = +E₀、m = m_e）与固定种子海面示意。
struct DiracSeaHoleModule: SimModule {

    let meta = ModuleMeta(
        id: "反粒子与狄拉克海_狄拉克海洞", title: "狄拉克海 · 洞 = 正电子",
        subtitle: "负能海填满 + 泡利缺位 → Q=+e、E=+m_e c² 的洞即正电子",
        category: .relativisticQFT, noteNumber: 39, tier: .realtime, difficulty: .basic,
        keywords: ["狄拉克海", "正电子", "反粒子", "负能态", "泡利原理", "洞",
                   "Dirac sea", "positron", "antiparticle"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "N", title: "海面示意点数", symbol: "N", unit: "",
                               range: 10...300, defaultValue: 90,
                               step: 10, decimalPlaces: 0)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .schematic(SchematicSpec(
                title: "狄拉克海：负能态填满，洞 = 正电子",
                xAxis: .init(label: "（示意横轴）"),
                yAxis: .init(label: "能量 E (MeV)")))
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let n = Int(input.slider("N"))
        let e0 = try DiracSeaMath.restEnergy(constants: constants)
        let eBound = 1.30

        let sea = DiracSeaMath.seaPoints(count: n, eBound: eBound)
        var seaMarkers: [SchematicMarker] = []
        seaMarkers.reserveCapacity(n)
        for p in sea {
            seaMarkers.append(.init(x: p.x, y: p.y, filled: true, colorIndex: 0))
        }

        return SimResult(
            charts: [
                .schematic(SchematicData(
                    spec: charts[0].schematicSpec!,
                    bands: [
                        .init(label: "负能连续谱（狄拉克海，已填满）",
                              xRange: -1.3...1.3, yRange: -eBound...0,
                              filled: true, colorIndex: 0),
                        .init(label: "正能连续谱（普通粒子）",
                              xRange: -1.3...1.3, yRange: 0...eBound,
                              filled: true, colorIndex: 3),
                    ],
                    markers: seaMarkers + [
                        .init(label: "缺失电子（洞）", x: 0, y: -e0,
                              filled: false, colorIndex: 2),
                        .init(label: "普通电子（Q = −e）", x: 0.6, y: 0.85,
                              filled: true, colorIndex: 0),
                        .init(label: "正电子（洞的显现）", x: 0.6, y: e0,
                              filled: true, colorIndex: 2),
                    ],
                    callouts: [
                        .init(text: "洞 = 正电子 e⁺\n(Q = +e, E = +\(String(format: "%.3f", e0)) MeV, m = m_e)",
                              anchor: Point(x: 0.95, y: e0),
                              arrowEnd: Point(x: 0.62, y: e0), colorIndex: 2),
                        .init(text: "缺失的负能电子\n(Q = −e, E = −\(String(format: "%.3f", e0)) MeV)",
                              anchor: Point(x: -0.95, y: -e0),
                              arrowEnd: Point(x: -0.05, y: -e0), colorIndex: 2),
                    ],
                    referenceLines: [
                        .init(label: "E = 0（海面费米能级）", axis: .y, value: 0, style: .threshold),
                    ])),
            ],
            summary: [
                .init(id: "e0", title: "洞能级 |E₀|",
                      value: String(format: "%.6f MeV", e0),
                      note: "m_e c²（CODATA \(constants.version.rawValue)）"),
                .init(id: "charge", title: "洞的电荷",
                      value: "+e", note: "满海总电荷 0 ⇒ 缺一个 −e 即 +e"),
                .init(id: "energy", title: "洞的能量（相对满海）",
                      value: "+0.511 MeV", note: "移走一个 E = −E₀ 电子使系统能量 +E₀"),
                .init(id: "mass", title: "洞的质量",
                      value: "m_e", note: "正电子与电子质量相同"),
                .init(id: "N", title: "海面示意点数",
                      value: "\(n)", note: "固定种子 20260901（脚本原值；SplitMix64 序列）"),
            ],
            theory: TheoryCard(
                title: "狄拉克海与洞（笔记 39 第六节 b）",
                formulas: [
                    "负能电子态（Q = −e，E = −E₀）按泡利原理全部填满",
                    "缺位（洞）相对满海携带：Q = +e，E = +E₀，m = m_e",
                    "⇒ 洞的行为 = 带正电、正能量、质量 m_e 的粒子 = 正电子",
                    "γ ≳ 2mc²：γ + 海中负能电子 → 电子 + 洞（正电子）",
                ],
                reading: "蓝色带是被负能电子（深点）填满的狄拉克海，黄色带是普通正能粒子。"
                    + "红色空心圈是海中缺失的一个负能电子——相对满海它等价于一个"
                    + "Q=+e、E=+m_ec² 的正电子（右下箭头）。这一「洞的符号翻转」"
                    + "正是狄拉克 1928 年预言反物质的机制。"))
    }
}