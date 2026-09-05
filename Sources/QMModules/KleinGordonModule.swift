import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/克莱因-戈登方程_色散关系.py）

/// 克莱因-戈登色散关系纯函数核。
/// E = ±√(p²c² + m²c⁴)；以 p·c（MeV）与 m₀c²（MeV）为单位即 E² = (pc)² + (m c²)²。
enum KGMath {

    /// 正能支 E₊(p) = √((pc)² + (mc²)²)（p·c 与 mc² 同单位）。
    static func energyPositive(_ pc: Double, m0: Double) -> Double {
        (pc * pc + m0 * m0).squareRoot()
    }

    /// 负能支 E₋(p) = −E₊(p)。
    static func energyNegative(_ pc: Double, m0: Double) -> Double {
        -energyPositive(pc, m0: m0)
    }

    /// 非相对论近似 E_NR = mc² + p²/2m（同一单位制）。
    static func energyNR(_ pc: Double, m0: Double) -> Double {
        m0 + pc * pc / (2.0 * m0)
    }

    /// 极端相对论极限下质量阈值占比 (mc²/E)²。
    static func massFraction(_ pc: Double, m0: Double) -> Double {
        let ratio = m0 / energyPositive(pc, m0: m0)
        return ratio * ratio
    }
}

// MARK: - 模块

/// 笔记 37《克莱因-戈登方程》：自由 KG 色散关系 E = ±√(p²c² + m²c⁴)
/// 的正负能两支与非相对论近似——p=0 处质量阈值 ±mc²。
struct KleinGordonDispersionModule: SimModule {

    let meta = ModuleMeta(
        id: "克莱因-戈登方程_色散关系", title: "克莱因-戈登 · 色散关系",
        subtitle: "E = ±√(p²c² + m²c⁴)——正负能两支与 p=0 质量阈值",
        category: .relativisticQFT, noteNumber: 37, tier: .realtime, difficulty: .basic,
        keywords: ["克莱因-戈登方程", "色散关系", "正能支", "负能支", "质量阈值",
                   "非相对论近似", "相对论量子力学", "Klein-Gordon", "dispersion"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "m0", title: "静止能量", symbol: "mc²", unit: "MeV",
                               range: 0.05...500, defaultValue: 0.51099895069,
                               scale: .log, decimalPlaces: 3)),
            .slider(SliderSpec(key: "pmax", title: "动量上限", symbol: "p·c", unit: "MeV",
                               range: 0.5...10, defaultValue: 2,
                               decimalPlaces: 1)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "KG 色散关系：正负能两支与非相对论近似",
                xAxis: .init(label: "p·c (MeV)"),
                yAxis: .init(label: "E (MeV)"),
                seriesNames: ["E₊ = +√((pc)² + (mc²)²)", "E₋ = −√((pc)² + (mc²)²)",
                              "非相对论 mc² + p²/2m"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let m0 = input.slider("m0")
        let pMax = input.slider("pmax")

        // 脚本同口径：p ∈ [−p_max, p_max] 2000 点
        let pc = Num.linspace(-pMax, pMax, count: 2000)
        let ePos = pc.map { KGMath.energyPositive($0, m0: m0) }
        let eNeg = pc.map { KGMath.energyNegative($0, m0: m0) }
        let eNR = pc.map { KGMath.energyNR($0, m0: m0) }

        // 关键数值（脚本 (1)-(4)）
        let pcHi = 10.0
        let eHi = KGMath.energyPositive(pcHi, m0: m0)
        let eHiNR = KGMath.energyNR(pcHi, m0: m0)

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: charts[0].lineSeriesSpec!,
                    series: [
                        .init(name: "E₊ = +√((pc)² + (mc²)²)",
                              points: Num.strided(pc, ePos, stride: 4)),
                        .init(name: "E₋ = −√((pc)² + (mc²)²)",
                              points: Num.strided(pc, eNeg, stride: 4), colorIndex: 1),
                        .init(name: "非相对论 mc² + p²/2m",
                              points: Num.strided(pc, eNR, stride: 4), colorIndex: 2),
                    ],
                    referenceLines: [
                        .init(label: "+mc² = \(String(format: "%.4f", m0)) MeV",
                              axis: .y, value: m0, style: .subtle),
                        .init(label: "−mc² = \(String(format: "%.4f", m0)) MeV",
                              axis: .y, value: -m0, style: .subtle),
                        .init(label: "p = 0", axis: .x, value: 0, style: .subtle),
                    ])),
            ],
            summary: [
                .init(id: "thresh", title: "p = 0 质量阈值",
                      value: String(format: "%.6f MeV", m0),
                      note: "E(0) = ±mc²——正负能支间隙 2mc²"),
                .init(id: "rel", title: "p·c = mc²（相对论性拐点）",
                      value: String(format: "%.6f MeV", KGMath.energyPositive(m0, m0: m0)),
                      note: "非相对论近似 \(String(format: "%.6f", KGMath.energyNR(m0, m0: m0))) MeV（开始失准）"),
                .init(id: "hi", title: "p·c = 10 MeV",
                      value: String(format: "%.6f MeV", eHi),
                      note: "非相对论近似 \(String(format: "%.6f", eHiNR)) MeV（高估）"),
                .init(id: "xr", title: "极端相对论质量占比",
                      value: String(format: "%.4e", KGMath.massFraction(pcHi, m0: m0)),
                      note: "(mc²/E)² → E ≈ |p|c"),
            ],
            theory: TheoryCard(
                title: "KG 色散关系（笔记 37）",
                formulas: [
                    "自由克莱因-戈登方程：E² = p²c² + m²c⁴",
                    "两支解：E₊ = +√(p²c² + m²c⁴)，E₋ = −√(p²c² + m²c⁴)",
                    "p = 0 处：E = ±mc²（质量阈值，正负能隙 2mc²）",
                    "非相对论近似：E ≈ mc² + p²/2m（p ≪ mc）；极端相对论 E ≈ |p|c（p ≫ mc）",
                ],
                reading: "蓝色正能支从 +mc² 单调抬升，红色负能支镜像下沉——两支在 p=0 处"
                    + "隔着 2mc² 的能隙，这是负能解（狄拉克海的伏笔）。绿色虚线是非相对论"
                    + "抛物线：只在 p ≪ mc 时贴合，p·c ≳ mc² 后显著高估能量。"))
    }
}