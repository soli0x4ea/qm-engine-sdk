import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/多电子原子_Zeff与电离能.py）

/// Slater 屏蔽规则 + NIST 第一电离能数据表。
enum MultiElectronAtomsMath {

    /// 壳层组态：(组标签, 电子数)。组标签首位数字编码 n，如 "2s2p" 表示 n=2 合并组。
    struct Element {
        let symbol: String
        let z: Int
        /// [(组标签, 电子数)]，Slater 顺序
        let config: [(group: String, count: Int)]
        /// 价电子所在组下标
        let target: Int
    }

    /// 前 19 号元素基态组态（与 Python elements 表逐行一致）。
    static let elements: [Element] = [
        .init(symbol: "H",  z: 1,  config: [("1s", 1)], target: 0),
        .init(symbol: "He", z: 2,  config: [("1s", 2)], target: 0),
        .init(symbol: "Li", z: 3,  config: [("1s", 2), ("2s2p", 1)], target: 1),
        .init(symbol: "Be", z: 4,  config: [("1s", 2), ("2s2p", 2)], target: 1),
        .init(symbol: "B",  z: 5,  config: [("1s", 2), ("2s2p", 3)], target: 1),
        .init(symbol: "C",  z: 6,  config: [("1s", 2), ("2s2p", 4)], target: 1),
        .init(symbol: "N",  z: 7,  config: [("1s", 2), ("2s2p", 5)], target: 1),
        .init(symbol: "O",  z: 8,  config: [("1s", 2), ("2s2p", 6)], target: 1),
        .init(symbol: "F",  z: 9,  config: [("1s", 2), ("2s2p", 7)], target: 1),
        .init(symbol: "Ne", z: 10, config: [("1s", 2), ("2s2p", 8)], target: 1),
        .init(symbol: "Na", z: 11, config: [("1s", 2), ("2s2p", 8), ("3s3p", 1)], target: 2),
        .init(symbol: "Mg", z: 12, config: [("1s", 2), ("2s2p", 8), ("3s3p", 2)], target: 2),
        .init(symbol: "Al", z: 13, config: [("1s", 2), ("2s2p", 8), ("3s3p", 3)], target: 2),
        .init(symbol: "Si", z: 14, config: [("1s", 2), ("2s2p", 8), ("3s3p", 4)], target: 2),
        .init(symbol: "P",  z: 15, config: [("1s", 2), ("2s2p", 8), ("3s3p", 5)], target: 2),
        .init(symbol: "S",  z: 16, config: [("1s", 2), ("2s2p", 8), ("3s3p", 6)], target: 2),
        .init(symbol: "Cl", z: 17, config: [("1s", 2), ("2s2p", 8), ("3s3p", 7)], target: 2),
        .init(symbol: "Ar", z: 18, config: [("1s", 2), ("2s2p", 8), ("3s3p", 8)], target: 2),
        .init(symbol: "K",  z: 19, config: [("1s", 2), ("2s2p", 8), ("3s3p", 8), ("4s4p", 1)], target: 3),
    ]

    /// Slater 屏蔽：同组其它电子 0.35（1s 为 0.30），n−1 层 0.85，n−2 及更深 1.00。
    static func slaterZeff(_ element: Element) -> Double {
        let config = element.config
        let nTarget = Int(String(config[element.target].group.prefix(1)))!
        var sigma = 0.0
        for (i, grp) in config.enumerated() {
            let nGrp = Int(String(grp.group.prefix(1)))!
            if i == element.target {
                let others = Double(grp.count - 1)
                sigma += others * (nGrp == 1 ? 0.30 : 0.35)
            } else if nGrp == nTarget - 1 {
                sigma += Double(grp.count) * 0.85
            } else if nGrp <= nTarget - 2 {
                sigma += Double(grp.count) * 1.00
            }
            // 更高 n 的组（右侧）按 Slater 规则贡献 0
        }
        return Double(element.z) - sigma
    }

    /// 全部元素的 Z_eff（按 Z 升序）。
    static var zeffList: [Double] { elements.map(slaterZeff) }

    /// 前 20 号元素第一电离能（eV）。
    /// 数据来源：NIST Atomic Spectra Database (ASD)，Kramida A. 等，
    /// NIST Standard Reference Database 78 —— 与源脚本同一数据表（单位 eV）。
    static let ionizationEnergies: [Double] = [
        13.598, 24.587, 5.392, 9.323, 8.298, 11.260, 14.534, 13.618,
        17.423, 21.565, 5.139, 7.646, 5.986, 8.152, 10.487, 10.360,
        12.968, 15.760, 4.341, 6.113,
    ]

    /// 稀有气体（峰）与碱金属（谷）的 Z 下标。
    static let nobleZ = [1, 9, 17]
    static let alkaliZ = [2, 10, 18]
}

// MARK: - 模块

/// 笔记 22《多电子原子与元素周期表》：Slater 屏蔽规则 Z_eff 趋势
/// + 前 20 号元素第一电离能周期性（NIST 实测）。实时档。
struct MultiElectronAtomsModule: SimModule {

    let meta = ModuleMeta(
        id: "多电子原子_Zeff与电离能", title: "多电子原子 · Zeff 与电离能",
        subtitle: "Slater 屏蔽规则 + 电离能周期性（NIST）",
        category: .manyBody, noteNumber: 22, tier: .realtime, difficulty: .basic,
        keywords: ["多电子", "屏蔽", "Slater", "有效核电荷", "Zeff", "Z_eff",
                   "电离能", "周期表", "壳层", "NIST"])

    var params: [ParamSpec] {
        [
            .discrete(DiscreteSpec(
                key: "element", title: "元素（价电子）",
                options: MultiElectronAtomsMath.elements.enumerated().map { i, e in
                    .init(id: e.symbol, title: e.symbol,
                          subtitle: "Z = \(e.z) · 组态 " +
                              e.config.map { "\($0.count)\($0.group)" }.joined())
                }, defaultOptionID: "Na")),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "Slater 有效核电荷 Z_eff（价电子）",
                xAxis: .init(label: "原子序数 Z"),
                yAxis: .init(label: "Z_eff"),
                seriesNames: ["Z_eff"])),
            .lineSeries(LineSeriesSpec(
                title: "第一电离能的周期性（NIST ASD）",
                xAxis: .init(label: "原子序数 Z"),
                yAxis: .init(label: "第一电离能 (eV)"),
                seriesNames: ["I₁(Z)", "稀有气体（峰）", "碱金属（谷）"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let symbol = input.discrete("element")
        let idx = MultiElectronAtomsMath.elements.firstIndex { $0.symbol == symbol } ?? 10
        let element = MultiElectronAtomsMath.elements[idx]
        let zeff = MultiElectronAtomsMath.slaterZeff(element)
        let zeffs = MultiElectronAtomsMath.zeffList
        let zSeq = MultiElectronAtomsMath.elements.map(\.z)
        let i1 = MultiElectronAtomsMath.ionizationEnergies
        let zFull = Array(1...20)

        let chart0 = LineSeriesData(
            spec: .init(xAxis: .init(label: "原子序数 Z"), yAxis: .init(label: "Z_eff"),
                        seriesNames: ["Z_eff"]),
            series: [.init(name: "Z_eff",
                           points: zip(zSeq, zeffs).map { Point(x: Double($0), y: $1) })],
            referenceLines: [ReferenceLine(
                id: "sel", label: String(format: "%@ Z_eff = %.2f", element.symbol, zeff),
                axis: .x, value: Double(element.z))])

        let noblePts = MultiElectronAtomsMath.nobleZ.map {
            Point(x: Double($0 + 1), y: i1[$0])
        }
        let alkaliPts = MultiElectronAtomsMath.alkaliZ.map {
            Point(x: Double($0 + 1), y: i1[$0])
        }
        let chart1 = LineSeriesData(
            spec: .init(xAxis: .init(label: "原子序数 Z"),
                        yAxis: .init(label: "第一电离能 (eV)"),
                        seriesNames: ["I₁(Z)", "稀有气体（峰）", "碱金属（谷）"]),
            series: [
                .init(name: "I₁(Z)",
                      points: zFull.map { Point(x: Double($0), y: i1[$0 - 1]) }),
                .init(name: "稀有气体（峰）", points: noblePts, colorIndex: 1),
                .init(name: "碱金属（谷）", points: alkaliPts, colorIndex: 2),
            ])

        let i1Sel = element.z <= 20 ? i1[element.z - 1] : 0
        return SimResult(
            charts: [.lineSeries(chart0), .lineSeries(chart1)],
            summary: [
                .init(id: "zeff", title: "Z_eff(\(element.symbol))",
                      value: String(format: "%.2f", zeff),
                      note: "Z=\(element.z) − σ (Slater)"),
                .init(id: "i1", title: "I₁(\(element.symbol))",
                      value: String(format: "%.3f eV", i1Sel),
                      note: "NIST ASD 实测值"),
                .init(id: "peaks", title: "峰/谷",
                      value: "He 24.59 · Ne 21.57 · Ar 15.76 / Li 5.39 · Na 5.14 · K 4.34 eV",
                      note: "稀有气体峰、碱金属谷"),
            ],
            theory: TheoryCard(
                title: "Slater 屏蔽与电离能周期性",
                formulas: [
                    "Z_eff = Z − σ：同组 0.35（1s 为 0.30），n−1 层 0.85，更深 1.00",
                    "壳层内 Z_eff 随 Z 线性上升 → 电离能上升",
                    "开新壳层（碱金属）σ 跳升 → Z_eff 骤降 → 电离能跌入谷",
                    "I₁ 数据：NIST Atomic Spectra Database (ASD)，SRD 78",
                ],
                reading: "Z_eff 曲线的锯齿就是周期表的壳层结构：每开新主壳层，"
                    + "价电子突然被完整内层屏蔽（Na 的 Z_eff 仅 2.2），电离能应声跌谷。"))
    }
}
