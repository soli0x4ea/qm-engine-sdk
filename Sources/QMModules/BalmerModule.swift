import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/原子光谱与玻尔模型_巴尔末系.py）

/// 玻尔模型巴尔末系的纯函数核（约化质量口径）。
enum BalmerMath {

    /// 约化质量因子 μ/m_e = m_p/(m_e + m_p)。
    static func reducedMassFactor(me: Double, mp: Double) -> Double {
        mp / (me + mp)
    }

    /// 氢里德伯常数 R_H = R_inf · μ/m_e。
    static func rydbergH(rInf: Double, me: Double, mp: Double) -> Double {
        rInf * reducedMassFactor(me: me, mp: mp)
    }

    /// 谱线波长（m）：1/λ = R_H(1/n1² − 1/n2²)。
    static func wavelength(n1: Int, n2: Int, rH: Double) -> Double {
        1.0 / (rH * (1.0 / Double(n1 * n1) - 1.0 / Double(n2 * n2)))
    }

    /// 能级 E_n = −E_ion/n²（eV，约化质量口径，E_ion = 13.598434 eV）。
    static func energyLevel(_ n: Int, eIon: Double = 13.598434) -> Double {
        -eIon / Double(n * n)
    }
}

// MARK: - 模块

/// 笔记 06《原子光谱与玻尔模型》：巴尔末系谱线 + 氢能级图（约化质量口径）。
/// 实时档：解析计算；多值对比组选择显示的谱线。
struct BalmerModule: SimModule {

    private static let lineDefs: [(id: String, label: String, n2: Int, obs: Double)] = [
        ("ha", "Hα (3→2)", 3, 656.469),
        ("hb", "Hβ (4→2)", 4, 486.273),
        ("hg", "Hγ (5→2)", 5, 434.173),
        ("hd", "Hδ (6→2)", 6, 410.293),
    ]

    let meta = ModuleMeta(
        id: "原子光谱与玻尔模型_巴尔末系", title: "原子光谱 · 巴尔末系",
        subtitle: "玻尔模型 Hα/Hβ/Hγ/Hδ 谱线与氢能级图",
        category: .oldQuantum, noteNumber: 6, tier: .realtime, difficulty: .basic,
        keywords: ["巴尔末", "Balmer", "玻尔", "Bohr", "里德伯", "Rydberg",
                   "氢原子", "光谱", "Hα", "能级", "跃迁"])

    var params: [ParamSpec] {
        [
            .multiCompare(MultiCompareSpec(
                key: "lines", title: "显示谱线",
                candidates: Self.lineDefs.map { .init(id: $0.id, label: $0.label, value: Double($0.n2)) },
                defaultSelectionIDs: Self.lineDefs.map(\.id), maxSelection: 4)),
            .constant(ConstantSpec(key: "R_inf", title: "R_∞", note: "里德伯常量（无限核质量）")),
            .constant(ConstantSpec(key: "m_e", title: "m_e", note: "电子质量")),
            .constant(ConstantSpec(key: "m_p", title: "m_p", note: "质子质量")),
        ]
    }

    var charts: [ChartSpec] {
        [
            .scatter(ScatterSpec(
                title: "巴尔末系（n₂→2，真空波长）",
                xAxis: .init(label: "λ (nm)"),
                yAxis: .init(label: "谱线"),
                seriesNames: ["计算值", "观测值"])),
            .levelDiagram(LevelDiagramSpec(
                title: "氢原子能级 E_n = −13.598 eV/n²",
                energyAxis: "E (eV)",
                showTransitions: true)),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let selected = input.multiCompare("lines")
        let rH = BalmerMath.rydbergH(rInf: try constants.value("R_inf"),
                                     me: try constants.value("m_e"),
                                     mp: try constants.value("m_p"))

        // 谱线散点：计算值 vs 真空观测标称值
        var calcPts: [Point] = [], obsPts: [Point] = []
        var summaries: [SummaryItem] = []
        for def in Self.lineDefs where selected.contains(def.id) {
            let lamNm = BalmerMath.wavelength(n1: 2, n2: def.n2, rH: rH) * 1e9
            calcPts.append(Point(x: lamNm, y: 1.0))
            obsPts.append(Point(x: def.obs, y: 0.97))
            let rel = abs(lamNm - def.obs) / def.obs * 100
            summaries.append(.init(id: def.id, title: def.label,
                                   value: String(format: "%.3f nm", lamNm),
                                   note: String(format: "观测 %.3f，偏差 %.3f%%", def.obs, rel)))
        }

        // 能级图：n=1..6 + 跃迁 2→1 / 3→2 / 4→2（与脚本 arrow 一致）
        let levels = (1...6).map {
            Level(label: "n=\($0)", energy: BalmerMath.energyLevel($0))
        }
        let transitions = [
            Transition(fromIndex: 1, toIndex: 0, label: "2→1"),
            Transition(fromIndex: 2, toIndex: 1, label: "3→2 (Hα)"),
            Transition(fromIndex: 3, toIndex: 1, label: "4→2 (Hβ)"),
        ]

        return SimResult(
            charts: [
                .scatter(ScatterData(
                    spec: ScatterSpec(
                        xAxis: .init(label: "λ (nm)"),
                        yAxis: .init(label: "谱线"),
                        seriesNames: ["计算值", "观测值"]),
                    series: [
                        .init(name: "计算值", points: calcPts, colorIndex: 0),
                        .init(name: "观测值", points: obsPts, colorIndex: 1),
                    ])),
                .levelDiagram(LevelDiagramData(
                    spec: LevelDiagramSpec(energyAxis: "E (eV)", showTransitions: true),
                    levels: levels, transitions: transitions)),
            ],
            summary: [.init(id: "rh", title: "R_H（约化质量）",
                            value: String(format: "%.3f m⁻¹", rH),
                            note: "R_∞·m_p/(m_e+m_p)")] + summaries,
            theory: TheoryCard(
                title: "巴尔末公式",
                formulas: [
                    "1/λ = R_H(1/2² − 1/n₂²)，n₂ = 3,4,5,6",
                    "R_H = R_∞·μ/m_e，μ = m_e·m_p/(m_e+m_p)",
                    "E_n = −13.598 eV/n²（约化质量口径）",
                    "Hα 656.5 nm：宇宙中最重要的发射线之一",
                ],
                reading: "切换谱线对照计算与观测值；能级图里 3→2 那一跳就是 Hα 红光的来源。"))
    }
}
