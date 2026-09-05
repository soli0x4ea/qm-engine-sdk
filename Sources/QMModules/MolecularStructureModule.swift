import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/分子结构_H2与苯.py）

/// H₂ Morse 势 + 苯 Hückel π 能级的纯函数核。
enum MolecularStructureMath {

    /// Morse 势 V(R) = D_e[(1 − e^(−a(R−R_e)))² − 1]，eV。E(R_e)=−D_e，E(∞)=0。
    static func morse(_ R: Double, de: Double, a: Double, re: Double) -> Double {
        let s = 1 - exp(-a * (R - re))
        return de * (s * s - 1)
    }

    /// 苯 Hückel π 能级：ε_k = α + 2βcos(2πk/6)（k = 0..5）。
    static func huckelLevels(beta: Double) -> [Double] {
        (0..<6).map { 2 * beta * cos(2 * .pi * Double($0) / 6) }
    }

    /// HOMO = α + β，LUMO = α − β（α 取 0 为参考）。
    static func homoLumoGap(beta: Double) -> (homo: Double, lumo: Double, gap: Double) {
        (beta, -beta, abs(-beta - beta))
    }
}

// MARK: - 模块

/// 笔记 23《分子结构与化学键》：H₂ Morse 势能井（实验参数）
/// + 苯 Hückel π 能级与 HOMO-LUMO 隙。实时档。
struct MolecularStructureModule: SimModule {

    let meta = ModuleMeta(
        id: "分子结构_H2与苯", title: "分子结构 · H₂ 与苯",
        subtitle: "Morse 势能井 + Hückel π 能级",
        category: .manyBody, noteNumber: 23, tier: .realtime, difficulty: .basic,
        keywords: ["分子", "化学键", "Morse", "莫尔斯", "H2", "氢分子", "苯",
                   "Hückel", "休克尔", "π电子", "HOMO", "LUMO", "离域"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "De", title: "解离能", symbol: "D_e", unit: "eV",
                               range: 1...8, defaultValue: 4.52, decimalPlaces: 2)),
            .slider(SliderSpec(key: "a", title: "宽度参数", symbol: "a", unit: "Å⁻¹",
                               range: 1...2.5, defaultValue: 1.94, decimalPlaces: 2)),
            .slider(SliderSpec(key: "Re", title: "平衡键长", symbol: "R_e", unit: "Å",
                               range: 0.5...1.2, defaultValue: 0.741, decimalPlaces: 3)),
            .slider(SliderSpec(key: "beta", title: "共振积分", symbol: "β", unit: "eV",
                               range: -4 ... -1, defaultValue: -2.5, decimalPlaces: 2)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "H₂ 键势能井（Morse，实验参数）",
                xAxis: .init(label: "核间距 R (Å)"),
                yAxis: .init(label: "势能 (eV)"),
                seriesNames: ["Morse V(R)"])),
            .levelDiagram(LevelDiagramSpec(
                title: "苯 Hückel π 能级（6 π 电子）",
                energyAxis: "能量 (eV, α=0)",
                showTransitions: false)),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let de = input.slider("De")
        let a = input.slider("a")
        let re = input.slider("Re")
        let beta = input.slider("beta")

        // 与 Python R = linspace(0.4, 4.0, 400) 一致（完整 400 点）
        let R = Num.linspace(0.4, 4.0, count: 400)
        let E = R.map { MolecularStructureMath.morse($0, de: de, a: a, re: re) }

        let chart0 = LineSeriesData(
            spec: .init(xAxis: .init(label: "核间距 R (Å)"),
                        yAxis: .init(label: "势能 (eV)"),
                        seriesNames: ["Morse V(R)"]),
            series: [.init(name: "Morse V(R)", points: zip(R, E).map { Point(x: $0, y: $1) })],
            referenceLines: [
                ReferenceLine(id: "re", label: String(format: "R_e = %.3f Å", re),
                              axis: .x, value: re),
                ReferenceLine(id: "de", label: String(format: "−D_e = %.2f eV", -de),
                              axis: .y, value: -de),
            ])

        let eps = MolecularStructureMath.huckelLevels(beta: beta)
        let occ = [true, true, false, false, false, true]  // k=0,1,5 占据（6 π 电子）
        let occNames = ["k=0", "k=1", "k=2", "k=3", "k=4", "k=5"]
        let levels = eps.enumerated()
            .sorted { $0.element < $1.element }
            .map { i, e in
                Level(label: "\(occNames[i]) " + (occ[i] ? "（占据）" : "（空）"),
                      energy: e)
            }
        let chart1 = LevelDiagramData(
            spec: .init(energyAxis: "能量 (eV, α=0)", showTransitions: false),
            levels: levels)

        let gap = MolecularStructureMath.homoLumoGap(beta: beta)
        return SimResult(
            charts: [.lineSeries(chart0), .levelDiagram(chart1)],
            summary: [
                .init(id: "well", title: "井深",
                      value: String(format: "E(R_e) = −%.2f eV", de),
                      note: "解离能 D_e（Herzberg 实验值 4.52 eV）"),
                .init(id: "homo", title: "HOMO",
                      value: String(format: "%+.2f eV", gap.homo),
                      note: "α + β"),
                .init(id: "lumo", title: "LUMO",
                      value: String(format: "%+.2f eV", gap.lumo),
                      note: "α − β"),
                .init(id: "gap", title: "π-π* 隙",
                      value: String(format: "%.2f eV", gap.gap),
                      note: "|LUMO − HOMO| = 2|β|"),
            ],
            theory: TheoryCard(
                title: "Morse 势与 Hückel 离域",
                formulas: [
                    "V(R) = D_e[(1 − e^(−a(R−R_e)))² − 1]：短程排斥 + 长程 −D_e 渐近",
                    "实验参数：D_e = 4.52 eV，a = 1.94 Å⁻¹，R_e = 0.741 Å",
                    "ε_k = α + 2βcos(2πk/6)：苯 6 环 6 π 电子的 Hückel 谱",
                    "占据：k=0（2e）+ k=1,5（4e）；HOMO = α+β，LUMO = α−β",
                ],
                reading: "Morse 井给出键长与解离能；苯的简并能级 2|β|=5 eV 隙对应"
                    + "紫外吸收——离域 π 电子正是苯环化学稳定性的来源。"))
    }
}
