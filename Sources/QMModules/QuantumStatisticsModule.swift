import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/量子统计_BEC与费米气体.py）

/// 理想玻色/费米气体的纯函数核。
enum QuantumStatisticsMath {

    /// BEC 临界温度 k_BT_c = (2πℏ²/m)(n/ζ(3/2))^(2/3)，ζ(3/2) = 2.612。
    static func becTc(n: Double, m: Double, hbar: Double, kB: Double) -> Double {
        let zeta32 = 2.612
        return 2 * .pi * hbar * hbar / m * pow(n / zeta32, 2.0 / 3) / kB
    }

    /// 凝聚分数 N₀/N = 1 − (T/T_c)^1.5（T < T_c），其余 0。
    static func condensateFraction(_ tRatio: Double) -> Double {
        tRatio < 1 ? 1 - pow(tRatio, 1.5) : 0
    }

    /// 费米能量 E_F = (ℏ²/2m)(3π²n)^(2/3)，J。
    static func fermiEnergy(_ n: Double, hbar: Double, me: Double) -> Double {
        hbar * hbar / (2 * me) * pow(3 * .pi * .pi * n, 2.0 / 3)
    }

    /// 零温简并压 P₀ = (2/5)nE_F（非相对论），Pa。
    static func degeneracyPressure(_ n: Double, eF: Double) -> Double {
        0.4 * n * eF
    }
}

// MARK: - 模块

/// 笔记 26《量子统计》：Rb-87 BEC 临界温度与凝聚分数 +
/// 费米简并压 P₀∝n^(5/3)（loglog）。实时档。
struct QuantumStatisticsModule: SimModule {

    let meta = ModuleMeta(
        id: "量子统计_BEC与费米气体", title: "量子统计 · BEC 与费米气体",
        subtitle: "凝聚分数 1−(T/T_c)^1.5 + 简并压 P₀∝n^(5/3)",
        category: .manyBody, noteNumber: 26, tier: .realtime, difficulty: .basic,
        keywords: ["量子统计", "玻色", "Bose", "BEC", "玻色-爱因斯坦凝聚", "临界温度",
                   "凝聚分数", "费米", "简并压", "白矮星", "Chandrasekhar"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "n_rb", title: "Rb-87 气体密度", symbol: "n",
                               unit: "m⁻³", range: 1e18...1e21, defaultValue: 1e19,
                               scale: .log, decimalPlaces: 2)),
            .slider(SliderSpec(key: "n_e", title: "电子密度", symbol: "n_e",
                               unit: "m⁻³", range: 1e27...1e36, defaultValue: 8.49e28,
                               scale: .log, decimalPlaces: 2)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "BEC 凝聚分数随温度",
                xAxis: .init(label: "T/T_c"),
                yAxis: .init(label: "凝聚分数 N₀/N"),
                seriesNames: ["N₀/N"])),
            .lineSeries(LineSeriesSpec(
                title: "零温费米简并压 vs 密度",
                xAxis: .init(label: "电子密度 n (m⁻³)", scale: .log),
                yAxis: .init(label: "简并压 P₀ (Pa)", scale: .log),
                seriesNames: ["P₀ ∝ n^(5/3)"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let nRb = input.slider("n_rb")
        let hbar = try constants.value("hbar")
        let kB = try constants.value("kB")
        let me = try constants.value("m_e")
        let eV = try constants.value("eV")
        let amu = try constants.value("u")

        // ---- 段 1：BEC（Rb-87，m = 86.909 u）----
        let mRb = 86.909 * amu
        let tc = QuantumStatisticsMath.becTc(n: nRb, m: mRb, hbar: hbar, kB: kB)
        let tRatio = Num.linspace(0, 1.3, count: 400)
        let frac = tRatio.map(QuantumStatisticsMath.condensateFraction)

        // ---- 段 2：费米简并压（n = logspace(27, 36, 400)）----
        let lnN = Num.linspace(27, 36, count: 400)
        let nGrid = lnN.map { pow(10, $0) }
        let p0Grid = nGrid.map {
            QuantumStatisticsMath.degeneracyPressure(
                $0, eF: QuantumStatisticsMath.fermiEnergy($0, hbar: hbar, me: me))
        }

        // 参考点：铜与白矮星（ρ = 1e9 kg/m³，μ_e = 2，非相对论近似）
        let nCu = 8.49e28
        let efCu = QuantumStatisticsMath.fermiEnergy(nCu, hbar: hbar, me: me)
        let p0Cu = QuantumStatisticsMath.degeneracyPressure(nCu, eF: efCu)
        let nWd = 1e9 / (2 * amu)
        let efWd = QuantumStatisticsMath.fermiEnergy(nWd, hbar: hbar, me: me)
        let p0Wd = QuantumStatisticsMath.degeneracyPressure(nWd, eF: efWd)

        let chart0 = LineSeriesData(
            spec: .init(xAxis: .init(label: "T/T_c"),
                        yAxis: .init(label: "凝聚分数 N₀/N"),
                        seriesNames: ["N₀/N"]),
            series: [.init(name: "N₀/N",
                           points: zip(tRatio, frac).map { Point(x: $0, y: $1) })],
            referenceLines: [ReferenceLine(label: "T/T_c = 1", axis: .x, value: 1)])

        let chart1 = LineSeriesData(
            spec: .init(xAxis: .init(label: "电子密度 n (m⁻³)", scale: .log),
                        yAxis: .init(label: "简并压 P₀ (Pa)", scale: .log),
                        seriesNames: ["P₀ ∝ n^(5/3)"]),
            series: [.init(name: "P₀ ∝ n^(5/3)",
                           points: zip(nGrid, p0Grid).map { Point(x: $0, y: $1) })],
            referenceLines: [
                ReferenceLine(id: "cu", label: "铜", axis: .x, value: nCu, style: .subtle),
                ReferenceLine(id: "wd", label: "白矮星（非相对论）",
                              axis: .x, value: nWd, style: .subtle),
            ])

        return SimResult(
            charts: [.lineSeries(chart0), .lineSeries(chart1)],
            summary: [
                .init(id: "tc", title: "T_c (Rb-87)",
                      value: String(format: "%.2e K（%.0f nK）", tc, tc * 1e9),
                      note: String(format: "n = %.2e m⁻³", nRb)),
                .init(id: "cu", title: "铜的简并压",
                      value: String(format: "E_F = %.3f eV，P₀ = %.3e Pa", efCu / eV, p0Cu),
                      note: "n = 8.49×10²⁸ m⁻³"),
                .init(id: "wd", title: "白矮星",
                      value: String(format: "E_F = %.0f keV，P₀ = %.3e Pa",
                                    efWd / eV / 1000, p0Wd),
                      note: "ρ = 10⁹ kg/m³，μ_e = 2（非相对论近似）"),
            ],
            theory: TheoryCard(
                title: "玻色凝聚与费米简并",
                formulas: [
                    "k_BT_c = (2πℏ²/m)·(n/ζ(3/2))^(2/3)，ζ(3/2) = 2.612",
                    "凝聚分数 N₀/N = 1 − (T/T_c)^(3/2)",
                    "E_F = (ℏ²/2m)(3π²n)^(2/3)；P₀ = (2/5)nE_F ∝ n^(5/3)",
                    "白矮星：ρ≈10⁹ kg/m³ → E_F≈164 keV——强简并（非相对论近似边缘）",
                ],
                reading: "同密度下 m 越大 T_c 越低（Rb-87 才 ~100 nK）；费米压不依赖温度，"
                    + "白矮星靠它抵抗引力——直到相对论性 P∝n^(4/3) 引发钱德拉塞卡极限。"))
    }
}
