import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/量子引力与全息原理_霍金温度.py / _黑洞熵.py）

/// 量子引力与全息原理纯函数核：霍金温度、黑洞熵（Bekenstein-Hawking）、普朗克单位。
enum QuantumGravityMath {
    /// 名义太阳质量（kg，脚本固定值，与 fixtures 一致）
    static let mSun: Double = 1.98847e30
    /// M87* 超大质量黑洞（≈ 6.5×10⁹ M☉）
    static let m87: Double = 6.5e9 * mSun

    static func planckMass(hbar: Double, c: Double, G: Double) -> Double {
        sqrt(hbar * c / G)
    }
    static func planckLength(hbar: Double, c: Double, G: Double) -> Double {
        sqrt(hbar * G / (c * c * c))
    }
    /// 霍金温度 T_H = ħ c³ / (8π G M k_B)
    static func hawkingTemperature(M: Double, hbar: Double, c: Double, G: Double, kB: Double) -> Double {
        hbar * c * c * c / (8 * .pi * G * M * kB)
    }
    /// 史瓦西半径 R_s = 2 G M / c²
    static func schwarzschildRadius(M: Double, c: Double, G: Double) -> Double {
        2 * G * M / (c * c)
    }
    /// Bekenstein-Hawking 熵（除以 k_B，无量纲计数）：S/k_B = A/(4 l_P²)，A = 4π R_s²
    static func bekensteinHawkingS(M: Double, hbar: Double, c: Double, G: Double) -> Double {
        let lP = planckLength(hbar: hbar, c: c, G: G)
        let Rs = schwarzschildRadius(M: M, c: c, G: G)
        return (4 * .pi * Rs * Rs) / (4 * lP * lP)
    }
}

// MARK: - 模块一：霍金温度

/// 笔记 41《量子引力与全息原理》：霍金温度
/// T_H = ħ c³ / (8π G M k_B)，log-log 下斜率为 −1（T_H ∝ M⁻¹）。
struct HawkingTemperatureModule: SimModule {

    let meta = ModuleMeta(
        id: "量子引力与全息原理_霍金温度", title: "霍金温度 · T_H ∝ M⁻¹",
        subtitle: "T_H = ħc³/(8πGMk_B)——黑洞越大越冷，log-log 斜率 −1",
        category: .relativisticQFT, noteNumber: 41, tier: .realtime, difficulty: .advanced,
        keywords: ["霍金温度", "Hawking", "黑洞", "事件视界", "量子引力", "Planck mass",
                   "black hole", "event horizon", "quantum gravity", "Hawking radiation"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "Mmax", title: "质量上限", symbol: "M_max", unit: "M☉",
                              range: 1e3...1e12, defaultValue: 1e10,
                              scale: .log, decimalPlaces: 0)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "霍金温度随黑洞质量变化（log-log）",
                xAxis: .init(label: "黑洞质量 M (kg)", scale: .log),
                yAxis: .init(label: "霍金温度 T_H (K)", scale: .log),
                seriesNames: ["T_H ∝ M^{-1}"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let hbar = try constants.value("hbar")
        let c = try constants.value("c")
        let G = try constants.value("G")
        let kB = try constants.value("kB")
        let mP = QuantumGravityMath.planckMass(hbar: hbar, c: c, G: G)
        let mSun = QuantumGravityMath.mSun
        let Mmax = input.slider("Mmax") * mSun

        let logM = Num.linspace(log10(mP), log10(Mmax), count: 400)
        let M = logM.map { pow(10.0, $0) }
        let T = M.map { QuantumGravityMath.hawkingTemperature(M: $0, hbar: hbar, c: c, G: G, kB: kB) }

        let refMarks: [(Double, String)] = [
            (mP, "Planck mass (2.18e-8 kg)"),
            (mSun, "1 M_sun"),
            (10 * mSun, "10 M_sun"),
            (QuantumGravityMath.m87, "M87*"),
        ]
        let refLines = refMarks.map {
            ReferenceLine(label: "\($0.1)\nT = \(QuantumGravityMath.hawkingTemperature(M: $0.0, hbar: hbar, c: c, G: G, kB: kB)) K",
                         axis: .x, value: $0.0, style: .subtle)
        }

        let tH_sun = QuantumGravityMath.hawkingTemperature(M: mSun, hbar: hbar, c: c, G: G, kB: kB)
        let tH_m87 = QuantumGravityMath.hawkingTemperature(M: QuantumGravityMath.m87, hbar: hbar, c: c, G: G, kB: kB)
        // W13h #41：log-log 直线自动域下拖 Mmax 画面不变——补端点标记 + 竖参考线
        let tHMax = QuantumGravityMath.hawkingTemperature(M: Mmax, hbar: hbar, c: c, G: G, kB: kB)

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: charts[0].lineSeriesSpec!,
                    series: [
                        .init(name: "T_H ∝ M^{-1}",
                              points: Num.strided(M, T, stride: 1)),
                    ],
                    referenceLines: refLines + [
                        ReferenceLine(
                            label: String(format: "M_max = %.2e M☉", input.slider("Mmax")),
                            axis: .x, value: Mmax, style: .threshold),
                    ],
                    pointMarkers: [
                        PointMarker(x: Mmax, y: tHMax,
                                    label: String(format: "M_max → T_H = %.3e K", tHMax),
                                    colorIndex: 1),
                    ])),
            ],
            summary: [
                .init(id: "mP", title: "普朗克质量 m_P",
                      value: String(format: "%.3e kg", mP),
                      note: "√(ħc/G)"),
                .init(id: "THsun", title: "太阳质量黑洞 T_H",
                      value: String(format: "%.3e K", tH_sun),
                      note: "≈ 6×10⁻⁸ K，远低于 CMB"),
                .init(id: "THm87", title: "M87* (6.5×10⁹ M☉) T_H",
                      value: String(format: "%.3e K", tH_m87),
                      note: "超大质量黑洞极冷"),
                .init(id: "slope", title: "log-log 斜率",
                      value: "−1", note: "T_H ∝ M⁻¹（ Hawking 标度律）"),
                .init(id: "mmax", title: "M_max 端点",
                      value: String(format: "%.2e M☉ → %.3e K",
                                    input.slider("Mmax"), tHMax),
                      note: "滑杆上限（图上大圆点，拖动即沿线滑动）"),
            ],
            theory: TheoryCard(
                title: "霍金温度（笔记 41 第二节、第六节）",
                formulas: [
                    "T_H = ħ c³ / (8π G M k_B)",
                    "log-log 斜率 = −1：黑洞质量越大，温度越低（T_H ∝ M⁻¹）",
                    "m_P = √(ħc/G) ≈ 2.18×10⁻⁸ kg；T_H(m_P) ≈ 平板温度上限",
                ],
                reading: "一条向下倾斜的直线：恒星级黑洞仅 ~10⁻⁸ K（比宇宙微波背景还冷，净吸收）；"
                    + "只有微观黑洞（接近普朗克质量）才灼烫。横轴标注了普朗克质量、太阳、M87* 等参考点，"
                    + "直观显示「越大越冷」的反直觉结论。"))
    }
}

// MARK: - 模块二：黑洞熵

/// 笔记 41《量子引力与全息原理》：Bekenstein-Hawking 黑洞熵
/// S_BH = A c³/(4Għ) = k_B A/(4 l_P²)，普朗克单位下 S/k_B = 4π(M/m_P)²；随质量单调增。
struct BlackHoleEntropyModule: SimModule {

    let meta = ModuleMeta(
        id: "量子引力与全息原理_黑洞熵", title: "黑洞熵 · S = A/4（普朗克单位）",
        subtitle: "S_BH = k_B A/(4 l_P²) = 4π(M/m_P)²——面积律与普朗克单位完全吻合",
        category: .relativisticQFT, noteNumber: 41, tier: .realtime, difficulty: .advanced,
        keywords: ["黑洞熵", "Bekenstein-Hawking", "面积律", "全息原理", "普朗克单位",
                   "black hole entropy", "holographic", "area law", "event horizon"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "Mmax", title: "质量上限", symbol: "M_max", unit: "M☉",
                              range: 1e1...1e12, defaultValue: 1e10,
                              scale: .log, decimalPlaces: 0)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "黑洞熵随质量变化（log-log）",
                xAxis: .init(label: "黑洞质量 M (kg)", scale: .log),
                yAxis: .init(label: "Bekenstein-Hawking 熵 S_BH / k_B", scale: .log),
                seriesNames: [
                    "S_BH/k_B = A/(4 l_P^2) ∝ M^2",
                    "Planck units: S_BH = 4π (M/m_P)^2",
                ])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let hbar = try constants.value("hbar")
        let c = try constants.value("c")
        let G = try constants.value("G")
        let mP = QuantumGravityMath.planckMass(hbar: hbar, c: c, G: G)
        let mSun = QuantumGravityMath.mSun
        let Mmax = input.slider("Mmax") * mSun

        // 主曲线：M 从 M_sun 起（恒星/超大质量黑洞；与脚本一致）
        let logM1 = Num.linspace(log10(mSun), log10(Mmax), count: 400)
        let M1 = logM1.map { pow(10.0, $0) }
        let S1 = M1.map { QuantumGravityMath.bekensteinHawkingS(M: $0, hbar: hbar, c: c, G: G) }

        // 普朗克单位对照曲线：M 从 m_P 起，理论上与主曲线完全重合
        let logM2 = Num.linspace(log10(mP), log10(Mmax), count: 400)
        let M2 = logM2.map { pow(10.0, $0) }
        let S2 = M2.map { 4 * .pi * pow($0 / mP, 2) }

        let refMarks: [(Double, String)] = [
            (mSun, "1 M_sun"),
            (QuantumGravityMath.m87, "M87*"),
        ]
        let refLines = refMarks.map {
            ReferenceLine(label: "\($0.1)\nS/k_B = \(QuantumGravityMath.bekensteinHawkingS(M: $0.0, hbar: hbar, c: c, G: G))",
                         axis: .x, value: $0.0, style: .subtle)
        }

        let sSun = QuantumGravityMath.bekensteinHawkingS(M: mSun, hbar: hbar, c: c, G: G)
        let sM87 = QuantumGravityMath.bekensteinHawkingS(M: QuantumGravityMath.m87, hbar: hbar, c: c, G: G)
        // W13h #41：同霍金温度——M_max 端点标记 + 竖参考线（拖滑杆沿线滑动）
        let sMax = QuantumGravityMath.bekensteinHawkingS(M: Mmax, hbar: hbar, c: c, G: G)

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: charts[0].lineSeriesSpec!,
                    series: [
                        .init(name: "S_BH/k_B = A/(4 l_P^2) ∝ M^2",
                              points: Num.strided(M1, S1, stride: 1)),
                        .init(name: "Planck units: S_BH = 4π (M/m_P)^2",
                              points: Num.strided(M2, S2, stride: 1), colorIndex: 1),
                    ],
                    referenceLines: refLines + [
                        ReferenceLine(
                            label: String(format: "M_max = %.2e M☉", input.slider("Mmax")),
                            axis: .x, value: Mmax, style: .threshold),
                    ],
                    pointMarkers: [
                        PointMarker(x: Mmax, y: sMax,
                                    label: String(format: "M_max → S/k_B = %.3e", sMax),
                                    colorIndex: 2),
                    ])),
            ],
            summary: [
                .init(id: "mP", title: "普朗克质量 m_P",
                      value: String(format: "%.3e kg", mP), note: "√(ħc/G)"),
                .init(id: "Ssun", title: "太阳质量黑洞 S/k_B",
                      value: String(format: "%.3e", sSun),
                      note: "≈ 10⁷⁷，远超可观测宇宙粒子数 ~10⁸⁰ 的量级框架"),
                .init(id: "Sm87", title: "M87* S/k_B",
                      value: String(format: "%.3e", sM87), note: "随 M² 暴涨"),
                .init(id: "law", title: "面积律",
                      value: "S = A/4", note: "S/k_B = A/(4 l_P²) = 4π(M/m_P)²（普朗克单位）"),
                .init(id: "mmax", title: "M_max 端点",
                      value: String(format: "%.2e M☉ → S/k_B = %.3e",
                                    input.slider("Mmax"), sMax),
                      note: "滑杆上限（图上大圆点，拖动即沿线滑动）"),
            ],
            theory: TheoryCard(
                title: "Bekenstein-Hawking 黑洞熵（笔记 41 第二节、第六节）",
                formulas: [
                    "S_BH = A c³/(4 G ħ) = k_B A/(4 l_P²)",
                    "A = 4π R_s²，R_s = 2 G M / c²",
                    "普朗克单位下：S/k_B = 4π (M/m_P)²（与面积律曲线完全重合）",
                ],
                reading: "熵正比于视界面积 A 而非体积——这是全息原理的雏形。绿线为普朗克单位表达式"
                    + " 4π(M/m_P)²，与蓝线（面积律）在重叠区间完全重合，验证两式等价。"
                    + "熵随质量平方暴涨：黑洞是熵密度最高的天体。"))
    }
}
