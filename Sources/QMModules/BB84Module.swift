import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/量子密钥分发_BB84密钥率.py）

/// BB84 渐近密钥率的纯函数核。
enum BB84Math {

    /// 二元香农熵 h(p) = −p log₂p − (1−p)log₂(1−p)（以 2 为底）。
    static func binaryEntropy(_ p: Double) -> Double {
        guard p > 0, p < 1 else { return 0 }
        return -(p * log2(p) + (1 - p) * log2(1 - p))
    }

    /// 单向后处理渐近密钥率 R(Q) = 1 − 2h(Q)（Shor-Preskill 2000）。
    static func keyRate(_ q: Double) -> Double {
        1 - 2 * binaryEntropy(q)
    }
}

// MARK: - 模块

/// 笔记 32《量子密钥分发》：BB84 单向渐近密钥率 R(Q) = 1 − 2h(Q)
/// 与 11% 安全阈值（区域填充标出安全区）。实时档。
struct BB84Module: SimModule {

    let meta = ModuleMeta(
        id: "量子密钥分发_BB84密钥率", title: "量子密钥分发 · BB84 密钥率",
        subtitle: "R(Q) = 1 − 2h(Q) 与 11% 安全阈值",
        category: .quantumInfo, noteNumber: 32, tier: .realtime, difficulty: .basic,
        keywords: ["量子密钥", "QKD", "BB84", "密钥率", "QBER", "误码率",
                   "Shor-Preskill", "香农熵", "安全阈值", "保密放大"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "Qmax", title: "横轴上限", symbol: "Q_max",
                               unit: "%", range: 15...40, defaultValue: 25,
                               decimalPlaces: 0)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "BB84 渐近密钥率 vs 误码率（单向后处理）",
                xAxis: .init(label: "QBER Q (%)"),
                yAxis: .init(label: "渐近密钥率 R (bit/pulse)"),
                seriesNames: ["R(Q) = 1 − 2h(Q)（单向）"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let qMax = input.slider("Qmax")

        // 与 Python Q = linspace(0, 0.25, 600) 同网格口径（按 Qmax 缩放）
        let qGrid = Num.linspace(0, qMax / 100, count: 600)
        let rGrid = qGrid.map(BB84Math.keyRate)

        // 安全区：Q ≤ 11%（固定阈值，与后处理方式绑定）
        let qOneWay = 0.11
        let rOneWay = BB84Math.keyRate(qOneWay)
        let safePts = qGrid.enumerated()
            .filter { $0.element <= qOneWay }
            .map { Point(x: $0.element * 100, y: max(rGrid[$0.offset], 0)) }

        let chart0 = LineSeriesData(
            spec: .init(xAxis: .init(label: "QBER Q (%)"),
                        yAxis: .init(label: "渐近密钥率 R (bit/pulse)"),
                        seriesNames: ["R(Q) = 1 − 2h(Q)（单向）"]),
            series: [.init(name: "R(Q) = 1 − 2h(Q)（单向）",
                           points: Num.strided(qGrid.map { $0 * 100 }, rGrid, stride: 2))],
            referenceLines: [
                ReferenceLine(id: "threshold", label: "单向阈值 Q = 11%",
                              axis: .x, value: 11),
                ReferenceLine(id: "zero", label: "R = 0", axis: .y, value: 0, style: .subtle),
            ],
            areaFills: [AreaFill(name: "安全区（R > 0）", points: safePts, baseline: 0)])

        let h11 = BB84Math.binaryEntropy(qOneWay)
        return SimResult(
            charts: [.lineSeries(chart0)],
            summary: [
                .init(id: "threshold", title: "单向安全阈值",
                      value: "Q* = 11.0%",
                      note: "h(Q*) = 0.49995 ≈ ½ 处 R=0"),
                .init(id: "r11", title: "R(11%)",
                      value: String(format: "%.4f", rOneWay),
                      note: "阈值处密钥率归零"),
                .init(id: "h", title: "h(0.11)",
                      value: String(format: "%.5f", h11),
                      note: "二元熵（bit）"),
                .init(id: "twoway", title: "双向后处理",
                      value: "阈值 ~12.6%", note: "优势蒸馏（不在本曲线内）"),
            ],
            theory: TheoryCard(
                title: "BB84 渐近密钥率（Shor-Preskill）",
                formulas: [
                    "R(Q) = 1 − 2h(Q)，h(Q) = −Q log₂Q − (1−Q)log₂(1−Q)",
                    "安全当且仅当 R > 0 ⟺ Q < 11.0%（单向后处理）",
                    "R(0) = 1：无误码时每脉冲提取 1 bit 密钥",
                    "双向后处理（优势蒸馏）阈值可升至 ~12.6%，速率公式不同",
                ],
                reading: "曲线在 Q=0 处满格 1 bit/pulse，随误码率二次塌缩"
                    + "（小 Q 时 h(Q)≈Q log₂(1/Q)）；填色区即「可提取密钥」的"
                    + "安全区——11% 右侧窃听者信息太多，保密放大吃光全部密钥。"))
    }
}
