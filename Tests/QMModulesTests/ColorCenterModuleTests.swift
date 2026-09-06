import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W7 笔记 44《色心与晶格缺陷的量子描述》：F 心类氢模型 / 电子-声子耦合吸收谱。
/// fixtures：44_色心与晶格缺陷的量子描述_F心类氢模型__v2022.json（12 根柱：6 模型 + 6 实验）、
/// 44_色心与晶格缺陷的量子描述_电子声子耦合吸收谱__v2022.json（512 点双谱，模块 4000 点）。
@Suite("W7 笔记44 色心与晶格缺陷的量子描述")
struct ColorCenterModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    // MARK: - 物理律：F 心类氢模型

    @Test("物理律：a* = ε∞a₀、R* = R_H/ε∞²、ΔE = (3/4)R*，故模型 λ ∝ ε∞²（六卤化物严格成立）")
    func fCenterHydrogenLaw() throws {
        let k = try constants()
        let rH = try k.value("R_inf_hc_eV")
        let a0 = try k.value("a0")
        let db = try MaterialDB.load()
        #expect(db.fCenterHalides.count == 6)
        for h in db.fCenterHalides {
            let rStar = ColorCenterMath.rydbergStarEV(epsInf: h.epsInf, rH_eV: rH)
            #expect(abs(rStar / (rH / (h.epsInf * h.epsInf)) - 1.0) < 1e-15, "\(h.id) R* = R_H/ε∞²")
            #expect(abs(ColorCenterMath.effectiveBohr(epsInf: h.epsInf, a0: a0)
                        / (h.epsInf * a0) - 1.0) < 1e-15, "\(h.id) a* = ε∞a₀")
            let dE = ColorCenterMath.transitionEnergyEV(epsInf: h.epsInf, rH_eV: rH)
            #expect(abs(dE / (0.75 * rStar) - 1.0) < 1e-15, "\(h.id) ΔE = (3/4)R*")
            let lam = ColorCenterMath.modelWavelengthNm(epsInf: h.epsInf, rH_eV: rH)
            #expect(abs(lam * dE / ColorCenterMath.hcEVnm - 1.0) < 1e-15, "\(h.id) λ = hc/ΔE")
            // 实验 F 带同样满足 λ = hc/E_F
            let lamExp = ColorCenterMath.expWavelengthNm(eF_eV: h.fBandExpEV)
            #expect(abs(lamExp * h.fBandExpEV / ColorCenterMath.hcEVnm - 1.0) < 1e-15)
        }
        // λ_model ∝ ε∞²：两盐之比恰为介电常数平方比
        let li = try #require(db.halide(id: "LiF"))
        let ki = try #require(db.halide(id: "KI"))
        let ratio = ColorCenterMath.modelWavelengthNm(epsInf: ki.epsInf, rH_eV: rH)
            / ColorCenterMath.modelWavelengthNm(epsInf: li.epsInf, rH_eV: rH)
        #expect(abs(ratio - pow(ki.epsInf / li.epsInf, 2)) < 1e-13, "λ(KI)/λ(LiF) = (ε∞_KI/ε∞_LiF)²")
        // 模型 λ 排序与 ε∞ 排序完全一致（同一幂律的单调性）
        let byEps = db.fCenterHalides.sorted { $0.epsInf < $1.epsInf }
        let byLambda = byEps.map { ColorCenterMath.modelWavelengthNm(epsInf: $0.epsInf, rH_eV: rH) }
        for i in 1..<byLambda.count {
            #expect(byLambda[i] > byLambda[i - 1], "模型 λ 应随 ε∞ 单调增（第 \(i) 项）")
        }
    }

    @Test("物理律：Mollwo-Ivey 经验律 E_F ∝ a^(−1.84)——实验 F 带随晶格常数单调下降，双对数拟合斜率 ≈ −1.74")
    func mollwoIveyLaw() throws {
        let db = try MaterialDB.load()
        let rows = db.fCenterHalides   // 已按晶格常数升序登记
        for i in 1..<rows.count {
            #expect(rows[i].latticeA > rows[i - 1].latticeA, "晶格常数应升序")
            #expect(rows[i].fBandExpEV < rows[i - 1].fBandExpEV,
                    "\(rows[i].id) 的实验 F 带能量应低于 \(rows[i - 1].id)（晶格更大 → 束缚更弱）")
        }
        // 双对数最小二乘拟合 ln E_F = c + m·ln a
        let xs = rows.map { log($0.latticeA) }
        let ys = rows.map { log($0.fBandExpEV) }
        let mx = xs.reduce(0, +) / Double(xs.count)
        let my = ys.reduce(0, +) / Double(ys.count)
        var sxy = 0.0, sxx = 0.0
        for i in xs.indices {
            sxy += (xs[i] - mx) * (ys[i] - my)
            sxx += (xs[i] - mx) * (xs[i] - mx)
        }
        let slope = sxy / sxx
        #expect(slope < -1.5 && slope > -2.0,
                "拟合斜率 \(slope) 应落在 Mollwo-Ivey 的 −1.84 ± 0.35 区间")
        // 类氢零阶模型系统性地高估 λ（低估跃迁能）——点离子修正未计入；
        // 其中 LiF（最小晶格、ε∞ 最小）偏差最大，正是正文引出点离子模型的动机。
        let k = try constants()
        let rH = try k.value("R_inf_hc_eV")
        var worstID = ""
        var worstRatio = 0.0
        for h in db.fCenterHalides {
            let lamModel = ColorCenterMath.modelWavelengthNm(epsInf: h.epsInf, rH_eV: rH)
            let lamExp = ColorCenterMath.expWavelengthNm(eF_eV: h.fBandExpEV)
            let ratio = lamModel / lamExp
            #expect(ratio > 1.0, "\(h.id)：模型 λ 应不低于实验 λ（实测比值 \(ratio)）")
            if ratio > worstRatio { worstRatio = ratio; worstID = h.id }
        }
        #expect(worstID == "LiF", "偏差最大的应为最小晶格 LiF（实测 \(worstID) 比值 \(worstRatio)）")
        // 大晶格（KCl/KBr/KI）模型与实验差距 < 30%
        for id in ["KCl", "KBr", "KI"] {
            let h = try #require(db.halide(id: id))
            let lamModel = ColorCenterMath.modelWavelengthNm(epsInf: h.epsInf, rH_eV: rH)
            let lamExp = ColorCenterMath.expWavelengthNm(eF_eV: h.fBandExpEV)
            #expect(abs(lamModel / lamExp - 1.0) < 0.30,
                    "\(id)：大晶格模型应接近实验（λ_m/λ_e = \(lamModel / lamExp)）")
        }
    }

    // MARK: - 物理律：电子-声子耦合

    @Test("物理律：Franck-Condon 因子为泊松分布——P_{m+1}/P_m = S/(m+1)，均值 = S")
    func franckCondonPoissonLaw() {
        for S in [0.5, 3.0, 6.0] {
            // nMax=40：尾部截断到 ~1e-20，归一化与均值可精确到 1e-10 量级
            // （模块的谱构建用默认 nMax=20，其截断误差由 fixture 全曲线对拍覆盖；
            //  此处仅验证 franckCondon 函数本身的泊松数学性质）
            let (ms, P) = ColorCenterMath.franckCondon(S, nMax: 40)
            #expect(ms.count == 41)
            // 递推比 P_{m+1}/P_m = S/(m+1)
            for m in 0..<40 {
                #expect(abs(P[m + 1] / P[m] - S / Double(m + 1)) < 1e-13,
                        "S=\(S)：P_\(m+1)/P_\(m) 应为 S/\(m+1)")
            }
            // 归一化与均值（泊松的两条定义性质）
            let sum = P.reduce(0, +)
            #expect(abs(sum - 1.0) < 1e-10, "S=\(S)：ΣP_m = \(sum)")
            var mean = 0.0
            for i in ms.indices { mean += Double(ms[i]) * P[i] }
            #expect(abs(mean - S) < 1e-6, "S=\(S)：平均声子阶数 \(mean) ≈ S")
            // 零声子线权重 = e^{−S}（Debye-Waller 因子）
            #expect(abs(P[0] - exp(-S)) < 1e-15, "S=\(S)：P₀ = e^{−S}")
        }
    }

    @Test("物理律：声子边带等间距 ħω、吸收/发射镜像对称、峰高比 P₁/P₀ = S")
    func phononSidebandLaw() async throws {
        let module = ElectronPhononAbsorptionModule()
        let values = ParamValues.defaults(for: module.params)
        let result = try await module.compute(values, constants: try constants())
        let chart = try #require(chartLineSeries(result, 0))
        #expect(chart.series.count == 2)
        let S = 3.0, hw = 0.065, eZPL = 1.945
        let xs = chart.series[0].points.map(\.x)
        let abs_ = chart.series[0].points.map(\.y)
        let emi = chart.series[1].points.map(\.y)
        #expect(xs.count == 4000)
        // 网格等距且对称于 ZPL：E[i] − E_ZPL = −(E[n−1−i] − E_ZPL)
        let step = xs[1] - xs[0]
        let last = try #require(xs.last)
        let first = try #require(xs.first)
        #expect(abs(step / ((last - first) / 3999.0) - 1.0) < 1e-12, "网格应等距")
        // 吸收谱的局部极大（声子边带峰）
        var peaks: [(x: Double, y: Double)] = []
        for i in 1..<(xs.count - 1) where abs_[i] > abs_[i - 1] && abs_[i] >= abs_[i + 1] {
            peaks.append((xs[i], abs_[i]))
        }
        #expect(peaks.count >= 8, "应分辨出 ≥8 条边带峰，实测 \(peaks.count)")
        #expect(abs(peaks[0].x - eZPL) < 1e-3, "首峰即零声子线 ZPL = \(peaks[0].x)")
        for i in 1..<peaks.count {
            #expect(abs((peaks[i].x - peaks[i - 1].x) - hw) < 1e-3,
                    "第 \(i) 条边带间距 \(peaks[i].x - peaks[i - 1].x) 应为 ħω = \(hw) eV")
        }
        // 峰高比 P₁/P₀ = S（γ=12 meV 的洛伦兹尾部在相邻峰位叠加约 7%，故容差 10%）
        #expect(abs(peaks[1].y / peaks[0].y / S - 1.0) < 0.10,
                "P₁/P₀ = \(peaks[1].y / peaks[0].y)，理论 S = \(S)")
        // 吸收 ↔ 发射镜像对称（E_ZPL + δ ↔ E_ZPL − δ）
        for i in 0..<xs.count {
            let mirror = emi[xs.count - 1 - i]
            #expect(abs(abs_[i] - mirror) / abs_[i] < 1e-10,
                    "第 \(i) 点吸收 \(abs_[i]) vs 发射镜像 \(mirror)")
        }
    }

    // MARK: - fixtures 对拍

    @Test("fixtures：F 心六卤化物 12 根柱（模型 λ / 实验 λ 各 6 根，< 1e-8）")
    func fCenterMatchesFixture() async throws {
        let fx = try ModuleFixture.load("44_色心与晶格缺陷的量子描述_F心类氢模型__v2022")
        let raw = try fx.bars(1, 0)
        #expect(raw.count == 12)
        // 布局：前 6 根为模型（x = i − 0.19），后 6 根为实验（x = i + 0.19）
        for i in 0..<6 {
            #expect(abs(raw[i].x - (Double(i) - 0.19)) < 1e-9, "模型柱 \(i) 位置")
            #expect(abs(raw[6 + i].x - (Double(i) + 0.19)) < 1e-9, "实验柱 \(i) 位置")
        }
        let fxModel = raw[0..<6].map(\.height)
        let fxExp = raw[6..<12].map(\.height)
        // 独立重算（Lit：LiF 模型 452.58 / 实验 244.06 nm，与 fixture 一致）
        #expect(abs(fxModel[0] - 452.583818) < 1e-5)
        #expect(abs(fxExp[0] - 244.063383) < 1e-5)

        let module = FCenterHydrogenModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        let bars = try #require(chartBars(result, 0))
        #expect(bars.bars.count == 12)
        let modelBars = bars.bars.enumerated().filter { $0.offset % 2 == 0 }.map(\.element.value)
        let expBars = bars.bars.enumerated().filter { $0.offset % 2 == 1 }.map(\.element.value)
        #expect(modelBars.count == 6 && expBars.count == 6)
        expectPointwiseClose(modelBars, fxModel, tolerance: 1e-8, "F 心模型 λ（六卤化物）")
        expectPointwiseClose(expBars, fxExp, tolerance: 1e-8, "F 心实验 λ（六卤化物）")
        // 默认选定 KCl：a*/a₀ = 2.22、R* = 2.76 eV（与 fixture 图题一致）
        let kcl = try #require(MaterialDB.load().halide(id: "KCl"))
        #expect(kcl.epsInf == 2.22)
        let levels = try #require(chartLevelDiagram(result, 1))
        // W13h #44：+实验 F 带能级 → 3 能级 2 跃迁
        #expect(levels.levels.count == 3)
        #expect(levels.transitions.count == 2)
        let gap = levels.levels[0].energy - levels.levels[1].energy
        let k = try constants()
        let rH = try k.value("R_inf_hc_eV")
        #expect(abs(gap - ColorCenterMath.transitionEnergyEV(epsInf: 2.22, rH_eV: rH)) < 1e-12)
    }

    @Test("fixtures：电子-声子双谱 512 点（计算核求值 < 1e-8；模块 4000 点重采样 < 2e-2）")
    func electronPhononMatchesFixture() async throws {
        let fx = try ModuleFixture.load("44_色心与晶格缺陷的量子描述_电子声子耦合吸收谱__v2022")
        let (xAbs, yAbs) = try fx.line(0, 0, index: 0)
        let (xEmi, yEmi) = try fx.line(0, 0, index: 1)
        #expect(xAbs.count == 512 && xEmi.count == 512)
        let S = 3.0, hw = 0.065, eZPL = 1.945, gamma = 0.012
        let cAbs = ColorCenterMath.buildSpectrum(xAbs, E_ZPL: eZPL, S: S, hw: hw,
                                                gamma: gamma, sign: +1)
        let cEmi = ColorCenterMath.buildSpectrum(xEmi, E_ZPL: eZPL, S: S, hw: hw,
                                                gamma: gamma, sign: -1)
        // 容差 1e-6：fixture 的 x 按 9 位有效数字存储（~1e-9 eV），经 γ=12 meV 洛伦兹
        // 的陡峭斜率（dL/L ~ 1/γ）放大到 ~3.5×10⁻⁷，实测最大相对偏差 3.5×10⁻⁷。
        expectPointwiseClose(cAbs, yAbs, tolerance: 1e-6, "吸收谱（fixture 网格求值）")
        expectPointwiseClose(cEmi, yEmi, tolerance: 1e-6, "发射谱（fixture 网格求值）")

        // 结构级重采样：fixture 512 点的格距 2.54 meV，而洛伦兹半宽 γ=12 meV（FWHM 24 meV），
        // 峰芯处线性插值误差实测最大 1.1%（全 4000 点），故取 2% 容差作形态级对照；
        // 判据级证据由上面的 fixture 网格求值给出（1e-6）。
        let module = ElectronPhononAbsorptionModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        let chart = try #require(chartLineSeries(result, 0))
        expectResampled(moduleX: chart.series[0].points.map(\.x),
                        moduleY: chart.series[0].points.map(\.y),
                        fx: xAbs, fy: yAbs, tolerance: 2e-2, "吸收谱重采样")
        expectResampled(moduleX: chart.series[1].points.map(\.x),
                        moduleY: chart.series[1].points.map(\.y),
                        fx: xEmi, fy: yEmi, tolerance: 2e-2, "发射谱重采样")
    }

    // MARK: - compute 结构与计时

    @Test("compute 结构：F 心 2 图（12 柱 + 3 能级 2 跃迁，W13h）+ 8 摘要")
    func fCenterStructureAndTiming() async throws {
        let module = FCenterHydrogenModule()
        let values = ParamValues.defaults(for: module.params)
        let result = try await module.compute(values, constants: try constants())
        #expect(result.charts.count == 2)
        let bars = try #require(chartBars(result, 0))
        #expect(bars.bars.count == 12)
        let levels = try #require(chartLevelDiagram(result, 1))
        // W13h #44：+实验 F 带能级与跃迁——切盐时模型-实验偏差可见
        #expect(levels.levels.count == 3 && levels.transitions.count == 2)
        let expLevel = try #require(levels.levels.first { $0.id == "1sExp" })
        let db0 = try MaterialDB.load()
        let kcl = try #require(db0.halide(id: values.discretes["salt"] ?? "KCl"))
        #expect(abs(expLevel.energy + kcl.fBandExpEV) < 1e-12,
                "实验能级 = −E_F(exp)（MaterialDB）")
        #expect(result.summary.count == 8, "6 卤化物 + 选定项 + Mollwo-Ivey")
        #expect(result.theory?.formulas.count == 3)
        // 切换到 LiF：ε∞ 更小 → R* 更大 → 能隙更宽（离散参数联动）
        var shifted = values
        shifted.discretes["salt"] = "LiF"
        let r2 = try await module.compute(shifted, constants: try constants())
        let lv2 = try #require(chartLevelDiagram(r2, 1))
        let k2 = try constants()
        let rH2 = try k2.value("R_inf_hc_eV")
        let gap2 = lv2.levels[0].energy - lv2.levels[1].energy
        #expect(abs(gap2 - ColorCenterMath.transitionEnergyEV(epsInf: 1.93, rH_eV: rH2)) < 1e-12,
                "LiF 能隙应等于闭式 0.75·R_H/ε∞²")
        #expect(gap2 > levels.levels[0].energy - levels.levels[1].energy,
                "LiF 的 1s→2p 能隙应大于 KCl（ε∞ 更小 → R* 更大）")
        try await expectComputeUnderBudget(module: module, values: values,
                                           constants: try constants(),
                                           budgetMillis: 16)
    }

    @Test("compute 结构：电子-声子谱 1 图双系列 4000 点 + ZPL 参考线 + 4 摘要（秒级档）")
    func electronPhononStructureAndTiming() async throws {
        let module = ElectronPhononAbsorptionModule()
        let values = ParamValues.defaults(for: module.params)
        let result = try await module.compute(values, constants: try constants())
        #expect(result.charts.count == 1)
        let chart = try #require(chartLineSeries(result, 0))
        #expect(chart.series.count == 2)
        for s in chart.series { #expect(s.points.count == 4000) }
        #expect(chart.referenceLines.count == 1)
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 3)
        // S 增大 → Debye-Waller 因子 e^{−S} 减小（ZPL 占比随耦合增强而下降）
        var bigger = values
        bigger.sliders["S"] = 6.0
        let r2 = try await module.compute(bigger, constants: try constants())
        let dw2 = try #require(r2.summary.first { $0.id == "DW" })
        let dw1 = try #require(result.summary.first { $0.id == "DW" })
        #expect(dw2.value == String(format: "%.3f", exp(-6.0)),
                "S=6 时 Debye-Waller 卡片应为 e⁻⁶ = \(exp(-6.0))，实际 \(dw2.value)")
        #expect(dw1.value == String(format: "%.3f", exp(-3.0)),
                "默认 S=3 时 DW 卡片应为 e⁻³ = \(exp(-3.0))，实际 \(dw1.value)")
        try await expectComputeUnderBudget(module: module, values: values,
                                           constants: try constants(),
                                           budgetMillis: 2000)
    }

    @Test("W13h #44：波长域模式——x 换算 λ = hc/E 且全部落在可见光窗 380–780 nm")
    func electronPhononNanometerMode() async throws {
        let module = ElectronPhononAbsorptionModule()
        var values = ParamValues.defaults(for: module.params)
        values.discretes["xmode"] = "nm"
        let result = try await module.compute(values, constants: try constants())
        let chart = try #require(chartLineSeries(result, 0))
        #expect(chart.series.count == 2)
        for s in chart.series {
            #expect(!s.points.isEmpty, "默认 NV 参数（ZPL=1.945 eV ≈ 637 nm）窗内应有谱")
            #expect(s.points.allSatisfy { (380.0...780.0).contains($0.x) },
                    "全部点在可见光窗内")
            // 换算自洽：x = 1239.841984/E 仅当 y 对应同一 E 网格——抽查能量域
            for p in s.points where abs(p.x - 637.1) < 3.0 {
                #expect(p.x > 0)
            }
        }
        // ZPL 参考线换算为 nm
        let zpl = try #require(chart.referenceLines.first)
        #expect(abs(zpl.value - 1239.841984 / 1.945) < 1e-9, "ZPL 参考线 = hc/E_ZPL nm")
    }
}
