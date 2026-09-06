import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/47_宝石学中的量子误用辨析_模型.py）

/// 宝石学量子误用辨析纯函数核：分离纯态 vs 纠缠约化最大混态的 Bloch 对照，
/// 红宝石 Cr³⁺ 简化能级（锐 R 发射 vs 蓝绿宽吸收）。3D 对照数据模型（SDK 交付物）：
/// 两态的 Bloch 矢量 + 半径 + 熵由 `dualBloch` 输出，SceneKit 渲染（剩余工作）直接消费。
enum GemologyMisuseMath {

    /// 分离纯态 |ψ(δ)⟩ = (|H⟩ + e^{iδ}|V⟩)/√2 的 Bloch 矢量。
    /// 密度矩阵 Pauli 迹路径（复用 BlochMath 态矢→r 核）：r = (cos δ, sin δ, 0)。
    static func separableBloch(delta: Double) -> (rx: Double, ry: Double, rz: Double) {
        let sv = BlochMath.stateVector(theta: .pi / 2, phi: delta)
        return BlochMath.blochVectorFromState(alphaRe: sv.alphaRe,
                                              betaRe: sv.betaRe, betaIm: sv.betaIm)
    }

    /// 纠缠 Bell 对单端约化态 ρ_A = I/2 的 Bloch 矢量：r = (0, 0, 0)。
    static func maximallyMixedBloch() -> (rx: Double, ry: Double, rz: Double) { (0, 0, 0) }

    // 红宝石 Cr³⁺ 能级（cm⁻¹；简化但以文献量级为准，脚本硬编码值）。
    static let e4A2 = 0.0          // 基态
    static let e2ER1 = 14402.0     // ≈ 694.3 nm（R1 锐发射线）
    static let e2ER2 = 14431.0     // ≈ 692.9 nm；R1/R2 零场分裂 ~29 cm⁻¹
    static let e4T2 = 18000.0      // 4A2→4T2 宽吸收 (~555 nm)——即八面体场 Δ_o
    static let e4T1 = 25000.0      // 4A2→4T1 宽吸收 (~400 nm)

    /// 2E 零场分裂（cm⁻¹）。
    static var zfs2E: Double { e2ER2 - e2ER1 }

    /// 波数 (cm⁻¹) → 真空波长 (nm)：λ = 10⁷ / ν̃。
    static func wavenumberToNm(_ cm: Double) -> Double { 1.0e7 / cm }

    /// 光子能量 (eV)：E = hc/λ，hc 取 CODATA 精确派生值（与 ColorMechanismMath 同一口径）。
    static let hcEVnm: Double = 1239.841984
    static func energyEV(wavelengthNm: Double) -> Double { hcEVnm / wavelengthNm }

    /// von Neumann 熵（bit）随 Bloch 半径 |r| 的二元熵 S = H₂((1+|r|)/2)。
    static func entropyVsRadius(_ radius: Double) -> Double {
        BlochMath.vonNeumannEntropyBits((radius, 0, 0))
    }
}

// MARK: - 模块：量子误用辨析模型（笔记 47 · 实时档）

/// 笔记 47《宝石学中的量子误用辨析》：双 Bloch 球对照（双折射后单光子纯态 |r|=1
/// vs 纠缠对单端约化最大混态 r=0 →「双折射不增纠缠」）+ 红宝石 Cr³⁺ 能级图
/// （锐 R 发射 vs 蓝绿宽吸收）。红宝石 Δ_o 消费 MaterialDB ruby 行做 round-trip。
struct GemologyMisuseModelModule: SimModule {

    let meta = ModuleMeta(
        id: "宝石学中的量子误用辨析_模型", title: "量子误用辨析 · Bloch 对照与红宝石能级",
        subtitle: "双折射不增纠缠（|r|=1 vs r=0，S=0 vs 1 bit）+ 红宝石 R 线 694.3 nm",
        category: .gemology, noteNumber: 47, tier: .realtime, difficulty: .advanced,
        keywords: ["量子误用", "纠缠", "混态", "纯态", "Bloch", "布洛赫", "红宝石", "Cr3+",
                   "R线", "694.3", "能级图", "双折射", "Bell", "约化态", "最大混态",
                   "ruby", "misuse", "entanglement"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "delta", title: "双折射相位差 δ", symbol: "δ", unit: "rad",
                              range: 0...(2 * .pi), defaultValue: .pi / 2,
                              scale: .linear, decimalPlaces: 3)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .bloch(BlochSpec(title: "Bloch 球对照：双折射纯态（赤道面随 δ 转动）")),
            .lineSeries(LineSeriesSpec(
                title: "von Neumann 熵 S 随 Bloch 半径 |r|（bit）",
                xAxis: .init(label: "Bloch 半径 |r|（0 = 球心最大混态，1 = 球面纯态）"),
                yAxis: .init(label: "纠缠熵 S (bit)"),
                seriesNames: ["S = H₂((1+|r|)/2)"])),
            .levelDiagram(LevelDiagramSpec(
                title: "红宝石 Cr³⁺ 能级：锐 R 发射 vs 蓝绿宽吸收",
                energyAxis: "能量 E (cm⁻¹)", showTransitions: true)),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let db = try MaterialDB.load()
        let delta = input.slider("delta")

        // 双 Bloch 对照（脚本图 1 的数据模型）
        let rSep = GemologyMisuseMath.separableBloch(delta: delta)
        let rEnt = GemologyMisuseMath.maximallyMixedBloch()
        let normSep = sqrt(rSep.rx * rSep.rx + rSep.ry * rSep.ry + rSep.rz * rSep.rz)
        let sSep = BlochMath.vonNeumannEntropyBits(rSep)
        let sEnt = BlochMath.vonNeumannEntropyBits(rEnt)

        // S(|r|) 曲线 + 球心/球面参考线
        let rGrid = Num.linspace(0.0, 1.0, count: 200)
        let sCurve = rGrid.map { GemologyMisuseMath.entropyVsRadius($0) }

        // 红宝石能级图（脚本图 2）
        let levels = [
            Level(label: "4A2 (ground)", energy: GemologyMisuseMath.e4A2),
            Level(label: "2E R1", energy: GemologyMisuseMath.e2ER1),
            Level(label: "2E R2", energy: GemologyMisuseMath.e2ER2),
            Level(label: "4T2", energy: GemologyMisuseMath.e4T2),
            Level(label: "4T1", energy: GemologyMisuseMath.e4T1),
        ]
        let r1Nm = GemologyMisuseMath.wavenumberToNm(GemologyMisuseMath.e2ER1)
        let r2Nm = GemologyMisuseMath.wavenumberToNm(GemologyMisuseMath.e2ER2)
        let t2Nm = GemologyMisuseMath.wavenumberToNm(GemologyMisuseMath.e4T2)
        let t1Nm = GemologyMisuseMath.wavenumberToNm(GemologyMisuseMath.e4T1)
        let transitions = [
            Transition(fromIndex: 1, toIndex: 0, label: String(format: "R1 发射 %.1f nm", r1Nm)),
            Transition(fromIndex: 2, toIndex: 0, label: String(format: "R2 发射 %.1f nm", r2Nm)),
            Transition(fromIndex: 0, toIndex: 3, label: String(format: "宽吸收 ~%.0f nm", t2Nm)),
            Transition(fromIndex: 0, toIndex: 4, label: String(format: "宽吸收 ~%.0f nm", t1Nm)),
        ]

        // MaterialDB round-trip：ruby 行 Δ_o 应等于 4T2 能级
        let ruby = db.mineral(id: "ruby")
        let deltaRuby = ruby?.deltaCm

        func fmt(_ v: Double) -> String { String(format: "%.5f", v) }

        return SimResult(
            charts: [
                .bloch(BlochData(
                    spec: charts[0].blochSpec!,
                    state: Point3D(x: rSep.rx, y: rSep.ry, z: rSep.rz),
                    stateLabel: "|ψ(δ)⟩",
                    trajectory: (0..<73).map { i in
                        let d = Double(i) / 72 * 2 * .pi
                        let r = GemologyMisuseMath.separableBloch(delta: d)
                        return Point3D(x: r.rx, y: r.ry, z: r.rz)
                    },
                    northLabel: "|R⟩", southLabel: "|L⟩")),
                .lineSeries(LineSeriesData(
                    spec: charts[1].lineSeriesSpec!,
                    series: [
                        .init(name: "S = H₂((1+|r|)/2)",
                              points: Num.strided(rGrid, sCurve, stride: 1)),
                    ],
                    referenceLines: [
                        ReferenceLine(label: "最大混态 |r|=0 → S=1", axis: .x, value: 0, style: .subtle),
                        ReferenceLine(label: "纯态 |r|=1 → S=0", axis: .x, value: 1, style: .subtle),
                    ],
                    pointMarkers: [
                        PointMarker(x: normSep, y: sSep,
                                    label: String(format: "δ=%.2f → |r|=%.2f", delta, normSep)),
                    ])),
                .levelDiagram(LevelDiagramData(spec: charts[2].levelDiagramSpec!,
                                               levels: levels, transitions: transitions)),
            ],
            summary: [
                .init(id: "sep", title: "双折射后单光子（可分离纯态）",
                      value: String(format: "r=(%.4f, %.4f, %.4f), |r|=%.4f",
                                    rSep.rx, rSep.ry, rSep.rz, normSep),
                      note: String(format: "δ=%.3f rad → r=(cos δ, sin δ, 0)；S=%.3f bit", delta, sSep)),
                .init(id: "mix", title: "纠缠 Bell 对单端约化态（最大混态）",
                      value: "r=(0, 0, 0), |r|=0",
                      note: String(format: "ρ_A = I/2；S=%.3f bit（1 bit）", sEnt)),
                .init(id: "radius", title: "径向对照：|r纯| − |r混| = 1",
                      value: String(format: "%.4f", normSep),
                      note: "「双折射不增纠缠」：幺正单比特演化把纯态转纯态，半径不变"),
                .init(id: "R1", title: "R1 锐发射线",
                      value: String(format: "%.1f nm (%.4f eV)", r1Nm,
                                    GemologyMisuseMath.energyEV(wavelengthNm: r1Nm)),
                      note: String(format: "2E R1 = %.0f cm⁻¹（脚本 694.3 nm 口径）",
                                   GemologyMisuseMath.e2ER1)),
                .init(id: "zfs", title: "2E 零场分裂 Δ",
                      value: String(format: "%.0f cm⁻¹", GemologyMisuseMath.zfs2E),
                      note: "R1/R2 双线间距（锐线 vs 宽带对照的特征）"),
                .init(id: "db", title: "MaterialDB ruby Δ_o round-trip",
                      value: deltaRuby.map { String(format: "%.0f cm⁻¹", $0) } ?? "缺失",
                      note: String(format: "= E(4T2) = %.0f cm⁻¹（宽吸收带即晶体场分裂）",
                                   GemologyMisuseMath.e4T2)),
            ],
            theory: TheoryCard(
                title: "量子误用辨析：纠缠、混态与红宝石发光（笔记 47）",
                formulas: [
                    "分离纯态 |ψ(δ)⟩ = (|H⟩+e^{iδ}|V⟩)/√2 → r = (cos δ, sin δ, 0)，|r|=1，S=0",
                    "纠缠对 |Φ⁺⟩ = (|HH⟩+|VV⟩)/√2 的单端约化态 ρ_A = I/2 → r=0，S=1 bit",
                    "双折射 = 单比特幺正演化：不改变约化密度矩阵谱 → 不能凭空增纠缠",
                    "红宝石：R 线锐发射 2E→4A2（694.3/692.9 nm，ΔZFS≈29 cm⁻¹）；蓝绿宽吸收 4A2→4T2/4T1",
                ],
                reading: "常见误用是把「偏振被双折射改变」说成「产生了纠缠」。Bloch 球上看一目了然："
                    + "单光子偏振态仍是球面上的纯态点（|r|=1），变的只是方向；而真正纠缠的 Bell 对，"
                    + "单端约化态是球心最大混态（r=0，S=1 bit）——纯与混的径向差恰为 1，这就是"
                    + "「双折射不增纠缠」的几何表述。红宝石的红色则来自 Cr³⁺ 能级结构：宽吸收带把"
                    + "蓝绿光泵入 4T2/4T1，无辐射弛豫到 2E，再沿锐 R 线（694.3 nm）发射——"
                    + "「宽吸收、锐发射」是激光与荧光宝石的能级图共性。"))
    }
}
