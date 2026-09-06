import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/44_色心与晶格缺陷的量子描述_F心类氢模型.py / _电子声子耦合吸收谱.py）

/// 色心与晶格缺陷纯函数核：F 心类氢有效质量近似、电子-声子吸收/发射谱（泊松边带）。
enum ColorCenterMath {
    /// h·c（eV·nm），CODATA 精确派生量（脚本固定 1239.841984）。
    static let hcEVnm: Double = 1239.841984

    /// 有效里德伯 R* = R_H / ε∞²（取 m*≈m_e）。
    static func rydbergStarEV(epsInf: Double, rH_eV: Double) -> Double {
        rH_eV / (epsInf * epsInf)
    }
    /// 有效玻尔半径 a* = ε∞ · a0。
    static func effectiveBohr(epsInf: Double, a0: Double) -> Double { epsInf * a0 }
    /// 1s→2p 跃迁能量 ΔE = (3/4) R*（即 0.75 R*）。
    static func transitionEnergyEV(epsInf: Double, rH_eV: Double) -> Double {
        0.75 * rydbergStarEV(epsInf: epsInf, rH_eV: rH_eV)
    }
    /// 类氢模型 F 带波长（nm）= hc / ΔE。
    static func modelWavelengthNm(epsInf: Double, rH_eV: Double) -> Double {
        hcEVnm / transitionEnergyEV(epsInf: epsInf, rH_eV: rH_eV)
    }
    /// 实验 F 带波长（nm）= hc / E_F(exp)。
    static func expWavelengthNm(eF_eV: Double) -> Double { hcEVnm / eF_eV }

    /// 归一化洛伦兹型（面积 = 1）：L(x) = (g/π) / ((x−x0)² + g²)。
    static func lorentzian(_ x: Double, _ x0: Double, _ g: Double) -> Double {
        (g / .pi) / ((x - x0) * (x - x0) + g * g)
    }

    /// Franck-Condon 因子 P_m = e^{−S} S^m / m!（泊松分布，T=0）。
    static func franckCondon(_ S: Double, nMax: Int) -> (ms: [Int], P: [Double]) {
        var ms: [Int] = [], P: [Double] = []
        let logBase = -S
        for m in 0...nMax {
            let logP = logBase + Double(m) * log(S) - lgamma(Double(m) + 1)
            ms.append(m); P.append(exp(logP))
        }
        return (ms, P)
    }

    /// 吸收(sign=+1)/发射(sign=−1) 谱：边带位于 E_ZPL ± m·ħω，每条洛伦兹线权重 P_m。
    static func buildSpectrum(_ E: [Double], E_ZPL: Double, S: Double, hw: Double,
                              gamma: Double, sign: Double, nMax: Int = 20) -> [Double] {
        let (ms, P) = franckCondon(S, nMax: nMax)
        return E.map { e in
            var s = 0.0
            for i in ms.indices {
                s += P[i] * lorentzian(e, E_ZPL + sign * Double(ms[i]) * hw, gamma)
            }
            return s
        }
    }
}

// MARK: - 模块一：F 心类氢模型

/// 笔记 44《色心与晶格缺陷的量子描述》：F 心（卤化物阴离子空位俘获电子）类氢有效质量近似
/// a* = ε∞ a0；R* = R_H/ε∞²；1s→2p 跃迁 ΔE = 0.75 R*。六卤化物对照（MaterialDB）。
struct FCenterHydrogenModule: SimModule {

    let meta = ModuleMeta(
        id: "色心与晶格缺陷的量子描述_F心类氢模型", title: "F 心 · 类氢模型六卤化物对照",
        subtitle: "a* = ε∞a₀；R* = R_H/ε∞²；ΔE = 0.75R*——六卤化物模型 vs 实验 F 带",
        category: .gemology, noteNumber: 44, tier: .realtime, difficulty: .advanced,
        keywords: ["F心", "色心", "类氢模型", "卤化物", "Mollwo-Ivey", "有效质量", "effective mass",
                   "color center", "F-center", "hydrogen-like", "alkali halide", "Mollwo-Ivey"])

    var params: [ParamSpec] {
        [
            .discrete(DiscreteSpec(
                key: "salt", title: "卤化物",
                options: (try? MaterialDB.load().fCenterHalides.map {
                    DiscreteOption(id: $0.id, title: $0.name)
                }) ?? [],
                defaultOptionID: "KCl")),
        ]
    }

    var charts: [ChartSpec] {
        [
            .bars(BarSpec(
                title: "F 心 1s→2p：类氢模型 λ vs 实验 F 带 λ（六卤化物）",
                xAxis: .init(label: "卤化物"),
                yAxis: .init(label: "F 带波长 λ (nm)"),
                seriesNames: ["Hydrogenic model λ", "Experimental F-band λ"])),
            .levelDiagram(LevelDiagramSpec(
                title: "F 心能级（选定卤化物，类氢近似）",
                energyAxis: "能量 E (eV)", showTransitions: true)),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let db = try MaterialDB.load()
        let rH = try constants.value("R_inf_hc_eV")
        let a0 = try constants.value("a0")

        var bars: [BarItem] = []
        var modelSummary: [SummaryItem] = []
        for h in db.fCenterHalides {
            let lamModel = ColorCenterMath.modelWavelengthNm(epsInf: h.epsInf, rH_eV: rH)
            let lamExp = ColorCenterMath.expWavelengthNm(eF_eV: h.fBandExpEV)
            bars.append(.init(label: h.name, value: lamModel, series: "Hydrogenic model λ"))
            bars.append(.init(label: h.name, value: lamExp, series: "Experimental F-band λ"))
            modelSummary.append(.init(
                id: "fc_\(h.id)", title: "\(h.name)：模型 λ / 实验 λ",
                value: String(format: "%.0f / %.0f nm", lamModel, lamExp),
                note: "ε∞=\(h.epsInf), E_F(exp)=\(h.fBandExpEV) eV"))
        }

        // 选定卤化物的能级图（W13h #44：加「实验 F 带」第三能级——
        // 原两能级图对 dE 自适应缩放，切卤化物画面恒不变；模型 vs 实验的
        // 相对偏差随盐切换可见，才是这页的物理）
        let selID = input.discrete("salt")
        let sel = db.halide(id: selID) ?? db.fCenterHalides[3]
        let rStar = ColorCenterMath.rydbergStarEV(epsInf: sel.epsInf, rH_eV: rH)
        let aStar = ColorCenterMath.effectiveBohr(epsInf: sel.epsInf, a0: a0)
        let dE = ColorCenterMath.transitionEnergyEV(epsInf: sel.epsInf, rH_eV: rH)
        let lam = ColorCenterMath.modelWavelengthNm(epsInf: sel.epsInf, rH_eV: rH)
        let dEexp = sel.fBandExpEV
        let lamExp = ColorCenterMath.expWavelengthNm(eF_eV: dEexp)
        let levels: [Level] = [
            .init(id: "2p", label: "2p (T1u)", energy: 0.0),
            .init(id: "1s", label: "1s 模型 (A1g)", energy: -dE),
            .init(id: "1sExp", label: "实验 F 带", energy: -dEexp),
        ]
        let transitions: [Transition] = [
            .init(fromIndex: 0, toIndex: 1,
                  label: "模型 ΔE = \(String(format: "%.2f", dE)) eV (λ≈\(Int(lam)) nm)"),
            .init(fromIndex: 0, toIndex: 2,
                  label: "实验 E_F = \(String(format: "%.2f", dEexp)) eV (λ≈\(Int(lamExp)) nm)"),
        ]

        return SimResult(
            charts: [
                .bars(BarData(spec: charts[0].barSpec!, bars: bars)),
                .levelDiagram(LevelDiagramData(spec: charts[1].levelDiagramSpec!,
                                              levels: levels, transitions: transitions)),
            ],
            summary: modelSummary + [
                .init(id: "sel", title: "选定：\(sel.name)",
                      value: String(format: "a* = %.2f a₀, R* = %.2f eV", aStar / a0, rStar),
                      note: "类氢零阶估计（大晶格如 KCl 接近实验）"),
                .init(id: "Mollwo", title: "Mollwo-Ivey 经验律",
                      value: "E_F ∝ a⁻¹·⁸⁴", note: "小晶格（LiF/NaCl）偏离大，引出点离子模型"),
            ],
            theory: TheoryCard(
                title: "F 心类氢有效质量近似（笔记 44 第六节）",
                formulas: [
                    "有效玻尔半径 a* = ε∞ a₀；有效里德伯 R* = R_H / ε∞²",
                    "1s→2p 跃迁 ΔE = (3/4) R* = 0.75 R* → λ = hc/ΔE",
                    "Mollwo-Ivey 经验律：E_F ∝ a⁻¹·⁸⁴（a 为晶格常数）",
                ],
                reading: "F 心是卤化物阴离子空位俘获一个电子形成的色心，可用类氢模型零阶估算其吸收带。"
                    + "模型波长随 ε∞² 增大而增大。对大晶格（如 KCl）模型接近实验，小晶格（LiF、NaCl）"
                    + "偏差显著——这正是正文引出点离子模型与 Mollwo-Ivey 经验律 E_F ∝ a⁻¹·⁸⁴ 的动机。"
                    + "下方柱状图为六卤化物模型 λ 与实验 F 带 λ 对照（数据来自 MaterialDB）。"))
    }
}

// MARK: - 模块二：电子-声子耦合吸收/发射谱

/// 笔记 44《色心与晶格缺陷的量子描述》：强电子-声子耦合下的吸收/发射谱
/// 零声子线 ZPL + 泊松权重的声子边带；Franck-Condon 因子 P_m = e^{−S} S^m/m!。
struct ElectronPhononAbsorptionModule: SimModule {

    let meta = ModuleMeta(
        id: "色心与晶格缺陷的量子描述_电子声子耦合吸收谱", title: "电子-声子谱 · ZPL + 泊松边带",
        subtitle: "吸收/发射谱 = Σ_m P_m·L(E_ZPL±mħω)；P_m 为泊松分布，边带间距 = ħω",
        category: .gemology, noteNumber: 44, tier: .seconds, difficulty: .advanced,
        keywords: ["电子声子耦合", "零声子线", "ZPL", "黄昆因子", "声子边带", "Debye-Waller",
                   "electron-phonon", "phonon sideband", "Huang-Rhys", "Franck-Condon", "NV center"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "S", title: "黄昆因子 S", symbol: "S", unit: "",
                              range: 0.5...10, defaultValue: 3.0,
                              scale: .linear, decimalPlaces: 2)),
            .slider(SliderSpec(key: "hw", title: "单声子能量 ħω", symbol: "ħω", unit: "eV",
                              range: 0.02...0.2, defaultValue: 0.065,
                              scale: .linear, decimalPlaces: 3)),
            .slider(SliderSpec(key: "EZPL", title: "零声子线能量", symbol: "E_ZPL", unit: "eV",
                              range: 1.0...3.0, defaultValue: 1.945,
                              scale: .linear, decimalPlaces: 3)),
            // W13h #44：X 轴可选波长域——NV 色心（钻石）谱在能量域一半落在红外，
            // 波长域下近紫外—可见区（300–700 nm，覆盖 NV ZPL 637 与激发侧特征）一目了然
            .discrete(DiscreteSpec(
                key: "xmode", title: "X 轴模式",
                options: [
                    .init(id: "ev", title: "能量域 E (eV)", subtitle: "脚本口径，ZPL ± 10ħω 全窗"),
                    .init(id: "nm", title: "波长域 λ (nm)", subtitle: "近紫外—可见区 300–700 nm 显示（少爷 W13j 口径）"),
                ],
                defaultOptionID: "ev")),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "电子-声子吸收/发射谱（T=0，ZPL + 泊松声子边带）",
                xAxis: .init(label: "光子能量 E (eV)"),
                yAxis: .init(label: "吸收/发射线型 (a.u.，log)", scale: .log),
                seriesNames: ["Absorption (T=0)", "Emission (T=0)"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let S = input.slider("S")
        let hw = input.slider("hw")
        let eZPL = input.slider("EZPL")
        let xmode = input.discrete("xmode")
        let gamma: Double = 0.012

        let nSide = 9
        let lo = eZPL - Double(nSide + 1) * hw
        let hi = eZPL + Double(nSide + 1) * hw
        let E = Num.linspace(lo, hi, count: 4000)
        let absSpec = ColorCenterMath.buildSpectrum(E, E_ZPL: eZPL, S: S, hw: hw, gamma: gamma, sign: +1)
        let emiSpec = ColorCenterMath.buildSpectrum(E, E_ZPL: eZPL, S: S, hw: hw, gamma: gamma, sign: -1)

        // W13h #44 + W13j #44：波长域模式——x 换算 λ = hc/E 并裁剪到 300–700 nm 窗
        let xAxisLabel: String
        var absPts = Num.strided(E, absSpec, stride: 1)
        var emiPts = Num.strided(E, emiSpec, stride: 1)
        var zplRef = eZPL
        if xmode == "nm" {
            func toNm(_ pts: [Point]) -> [Point] {
                pts.compactMap { p in
                    p.x > 0 ? Point(x: ColorCenterMath.hcEVnm / p.x, y: p.y) : nil
                }
                .filter { (300.0...700.0).contains($0.x) }
            }
            absPts = toNm(absPts)
            emiPts = toNm(emiPts)
            zplRef = ColorCenterMath.hcEVnm / eZPL
            xAxisLabel = "波长 λ (nm)（近紫外—可见 300–700）"
        } else {
            xAxisLabel = "光子能量 E (eV)"
        }

        let wTL0 = exp(-S)                 // 零温 Debye-Waller 因子 = P_0
        let zMask = E.indices.filter { abs(E[$0] - eZPL) < gamma * 2.5 }
        let areaTotal = Num.trapezoid(E, absSpec)
        let areaZPL = Num.trapezoid(zMask.map { E[$0] }, zMask.map { absSpec[$0] })

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: .init(xAxis: .init(label: xAxisLabel),
                                yAxis: .init(label: "吸收/发射线型 (a.u.，log)", scale: .log),
                                seriesNames: ["Absorption (T=0)", "Emission (T=0)"]),
                    series: [
                        .init(name: "Absorption (T=0)", points: absPts),
                        .init(name: "Emission (T=0)", points: emiPts, colorIndex: 1),
                    ],
                    referenceLines: [
                        ReferenceLine(label: xmode == "nm"
                            ? "ZPL = \(String(format: "%.0f", zplRef)) nm"
                            : "ZPL = \(String(format: "%.3f", zplRef)) eV",
                                     axis: .x, value: zplRef, style: .threshold),
                    ])),
            ],
            summary: [
                .init(id: "S", title: "黄昆因子 S", value: String(format: "%.2f", S),
                      note: "强耦合（NV 中心典型 ~3）"),
                .init(id: "DW", title: "零温 Debye-Waller 因子 e^{−S}",
                      value: String(format: "%.3f", wTL0),
                      note: "= P₀，ZPL 占总强度比例"),
                .init(id: "DWint", title: "数值积分 ZPL 占比",
                      value: String(format: "%.3f", areaZPL / max(areaTotal, 1e-300)),
                      note: "应与 e^{−S} 一致"),
                .init(id: "hw", title: "边带间距 ħω", value: String(format: "%.0f meV", hw * 1000),
                      note: "相邻声子边带能量差"),
            ],
            theory: TheoryCard(
                title: "电子-声子耦合吸收谱（笔记 44 第二节、第六节）",
                formulas: [
                    "H_ep 黄昆因子 S；Franck-Condon 因子 P_m = e^{−S} S^m/m!（泊松）",
                    "吸收谱：边带位于 E_ZPL + m·ħω；发射谱：E_ZPL − m·ħω",
                    "零温 Debye-Waller 因子 W = e^{−S} = P₀（ZPL 相对强度）",
                ],
                reading: "强电子-声子耦合下，电子跃迁伴随晶格振动量子（声子）的发射/吸收，"
                    + "形成以零声子线 ZPL 为中心、间距恰为单声子能量 ħω 的泊松权重边带。"
                    + "黄昆因子 S 越大，振动弛豫越强、ZPL 越弱（Debye-Waller 因子 e^{−S}）。"
                    + "以 NV 中心（S≈3, ħω≈65 meV）为例，蓝线吸收、红线发射对称分布于 ZPL 两侧。"))
    }
}
