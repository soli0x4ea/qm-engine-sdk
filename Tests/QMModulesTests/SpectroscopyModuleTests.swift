import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W7 笔记 45《光谱学仪器的量子基础》：爱因斯坦 A 系数与荧光寿命 / 拉曼 Stokes–anti-Stokes 谱。
/// fixtures：45_光谱学仪器的量子基础_荧光与寿命__v2022.json（衰减双线 400 点，同网格；
/// 双轴 A(λ)/τ(λ) 图未导出数据，仅结构核对）、45_光谱学仪器的量子基础_拉曼谱__v2022.json
///（512 点，模块 60001 点 → 重采样对拍）。
@Suite("W7 笔记45 光谱学仪器的量子基础")
struct SpectroscopyModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    /// 脚本固定的 H 2p→1s 标定：d = 0.7449 e·a₀，λ = 121.567 nm（Lyman-α）
    private let dH_factor = 0.7449
    private let lamH = 121.567e-9

    private func hydrogenRates() throws -> (A: Double, tau: Double) {
        let k = try constants()
        let dH = dH_factor * (try k.value("e")) * (try k.value("a0"))
        let A = SpectroscopyMath.einsteinA(wavelength_m: lamH, d_Cm: dH,
                                          h: try k.value("h"), eps0: try k.value("epsilon_0"),
                                          hbar: try k.value("hbar"), c: try k.value("c"))
        return (A, 1.0 / A)
    }

    // MARK: - 物理律：爱因斯坦 A 系数与寿命

    @Test("物理律：A ∝ λ⁻³（ω³ 标度）且 τ = 1/A——双轴两系列逐点互为倒数")
    func einsteinALaw() async throws {
        let module = FluorescenceLifetimeModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        let chart = try #require(chartDualAxis(result, 0))
        #expect(chart.primary.count == 1 && chart.secondary.count == 1)
        let aSeries = chart.primary[0].points
        let tSeries = chart.secondary[0].points
        #expect(aSeries.count == 300 && tSeries.count == 300)
        // 共享 x 网格（发射波长 nm）
        for i in aSeries.indices {
            #expect(aSeries[i].x == tSeries[i].x, "第 \(i) 点 x 网格应一致")
            // τ(ns) = 1/A × 1e9 ⟺ A·τ = 1
            #expect(abs(aSeries[i].y * tSeries[i].y * 1e-9 - 1.0) < 1e-12,
                    "第 \(i) 点：A·τ = \(aSeries[i].y * tSeries[i].y * 1e-9)")
        }
        // ω³ 标度：λ 加倍 → A 变为 1/8
        for i in 1..<aSeries.count {
            let ratio = aSeries[i].y / aSeries[i - 1].y
            let expected = pow(aSeries[i - 1].x / aSeries[i].x, 3.0)
            #expect(abs(ratio / expected - 1.0) < 1e-12,
                    "第 \(i) 点：A 比 \(ratio) 应为 (λ⁻³) \(expected)")
        }
        // 短波辐射更快：A(300 nm) > A(800 nm)，τ 反之
        #expect(aSeries[0].y > aSeries[299].y)
        #expect(tSeries[0].y < tSeries[299].y)
    }

    @Test("物理律：H 2p→1s 标定——A ≈ 6.26×10⁸ s⁻¹、τ ≈ 1.6 ns（Lyman-α 公认寿命）")
    func lymanAlphaCalibration() async throws {
        let (A, tau) = try hydrogenRates()
        #expect(A > 6.0e8 && A < 6.5e8, "A_H = \(A) s⁻¹（文献 ~6.26×10⁸）")
        #expect(tau * 1e9 > 1.5 && tau * 1e9 < 1.7, "τ_H = \(tau * 1e9) ns（文献 ~1.6 ns）")
        // A ∝ |d|²：偶极矩加倍 → A 四倍
        let k = try constants()
        let d0 = 1.0 * SpectroscopyMath.debyeCm
        let a1 = SpectroscopyMath.einsteinA(wavelength_m: 500e-9, d_Cm: d0,
                                           h: try k.value("h"), eps0: try k.value("epsilon_0"),
                                           hbar: try k.value("hbar"), c: try k.value("c"))
        let a2 = SpectroscopyMath.einsteinA(wavelength_m: 500e-9, d_Cm: 2 * d0,
                                           h: try k.value("h"), eps0: try k.value("epsilon_0"),
                                           hbar: try k.value("hbar"), c: try k.value("c"))
        #expect(abs(a2 / a1 - 4.0) < 1e-15, "A ∝ |d|²")
        // 模块摘要的 τ_fluor 应等于本标定（ns 级允许自旋禁戒跃迁的口径差）
        let module = FluorescenceLifetimeModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        let card = try #require(result.summary.first { $0.id == "tauFluor" })
        #expect(card.value == String(format: "%.2f ns", tau * 1e9),
                "摘要 τ_fluor = \(card.value)，标定值 \(tau * 1e9) ns")
    }

    @Test("物理律：荧光/磷光双指数衰减——exp(−t/τ) 恒等、单调、寿命差 ~6 个数量级")
    func decayLaw() async throws {
        let (_, tauH) = try hydrogenRates()
        let module = FluorescenceLifetimeModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        let chart = try #require(chartLineSeries(result, 1))
        #expect(chart.series.count == 2)
        let t = chart.series[0].points.map(\.x)
        let iF = chart.series[0].points.map(\.y)
        let iP = chart.series[1].points.map(\.y)
        #expect(t.count == 400)
        let tauPhos = 3.0e-3
        for i in t.indices {
            #expect(abs(iP[i] - exp(-t[i] / tauPhos)) < 1e-12, "磷光第 \(i) 点")
            #expect(iP[i] > 0 && iP[i] <= 1.0)
            if iF[i] > 1e-12 {
                #expect(abs(iF[i] - exp(-t[i] / tauH)) / iF[i] < 1e-9, "荧光第 \(i) 点")
            }
            if i > 0 {
                #expect(iF[i] <= iF[i - 1] && iP[i] <= iP[i - 1], "衰减应单调不增")
            }
        }
        // 寿命比 ~1.9×10⁶（磷光自旋禁戒，慢六个数量级）
        let ratio = tauPhos / tauH
        #expect(ratio > 1e6 && ratio < 1e7, "τ_phos/τ_fluor = \(ratio)")
        let ratioCard = try #require(result.summary.first { $0.id == "ratio" })
        #expect(ratioCard.value == String(format: "%.1e", ratio))
    }

    // MARK: - 物理律：拉曼散射

    @Test("物理律：Stokes/anti-Stokes 峰位 = ν₀ ∓ ν_v，强度比 = (ν_AS/ν_S)⁴·e^{−ħω/kT}")
    func ramanBoltzmannLaw() async throws {
        let k = try constants()
        let c = try k.value("c"), hbar = try k.value("hbar"), kB = try k.value("kB")
        let c_cm = c * 100.0
        let T = 300.0
        let db = try MaterialDB.load()
        let modes = db.ramanModes
        #expect(modes.count == 3)
        let module = RamanSpectrumModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        let chart = try #require(chartLineSeries(result, 0))
        #expect(chart.series[0].points.count == 60001)
        let xs = chart.series[0].points.map(\.x)
        let ys = chart.series[0].points.map(\.y)
        let firstX = try #require(xs.first)
        let lastX = try #require(xs.last)
        let nu0 = (firstX + lastX) / 2.0
        #expect(abs(nu0 - 1.0e7 / 785.0) < 1e-6, "ν₀ = \(nu0) cm⁻¹（785 nm 激发）")

        // 局部极大 = 拉曼带（3 个模式 × Stokes/anti-Stokes = 6 条）
        var peaks: [(x: Double, y: Double)] = []
        for i in 1..<(xs.count - 1) where ys[i] > ys[i - 1] && ys[i] >= ys[i + 1] {
            peaks.append((xs[i], ys[i]))
        }
        #expect(peaks.count == 6, "应分辨 6 条拉曼带，实测 \(peaks.count)")
        let stokes = peaks.filter { $0.x < nu0 }.sorted { $0.x > $1.x }   // ν₀−ν_v，ν_v 升序
        let anti = peaks.filter { $0.x > nu0 }.sorted { $0.x < $1.x }
        #expect(stokes.count == 3 && anti.count == 3)
        for i in 0..<3 {
            let nuV = modes.sorted { $0.nuVcm < $1.nuVcm }[i].nuVcm
            #expect(abs((nu0 - stokes[i].x) - nuV) < 0.2,
                    "第 \(i) 条 Stokes 位移 \(nu0 - stokes[i].x) 应为 ν_v = \(nuV)")
            #expect(abs((anti[i].x - nu0) - nuV) < 0.2,
                    "第 \(i) 条 anti-Stokes 位移 \(anti[i].x - nu0) 应为 ν_v = \(nuV)")
            // 强度比：anti-Stokes/Stokes = (ν_AS/ν_S)⁴·n_v/(n_v+1)。
            // 弱 anti-Stokes 峰（尤其 ν_v=2900）坐在最强模式的洛伦兹拖尾上，直接取模块谱
            // 局部极大值会把拖尾算进去；这里减去「其它模式拖尾 + 本模 Stokes 拖尾」，
            // 还原单模纯振幅再比，物理律判据才干净。
            let omegaV = 2.0 * .pi * c_cm * nuV
            let nv = 1.0 / (exp(hbar * omegaV / (kB * T)) - 1.0)
            let nuS = nu0 - nuV, nuAS = nu0 + nuV
            let aSi = stokes[i].y                          // ≈ 本模 Stokes 纯振幅（污染 < 1e-5）
            let fwhm_i = try #require(modes.first { $0.nuVcm == nuV }).fwhmCm
            let h_i = fwhm_i / 2.0
            // 本模 Stokes 拖尾（在 anti 峰位处）
            var tail = aSi * (h_i * h_i) / ((anti[i].x - nuS) * (anti[i].x - nuS) + h_i * h_i)
            // 其它模式在 anti 峰位处的 Stokes + anti-Stokes 拖尾
            for mj in modes where mj.nuVcm != nuV {
                let oj = 2.0 * .pi * c_cm * mj.nuVcm
                let nvj = 1.0 / (exp(hbar * oj / (kB * T)) - 1.0)
                let nS = nu0 - mj.nuVcm, nAS = nu0 + mj.nuVcm
                let aSj = pow(nS, 4) * (nvj + 1.0), aASj = pow(nAS, 4) * nvj
                let hj = mj.fwhmCm / 2.0
                tail += aSj * (hj * hj) / ((anti[i].x - nS) * (anti[i].x - nS) + hj * hj)
                      + aASj * (hj * hj) / ((anti[i].x - nAS) * (anti[i].x - nAS) + hj * hj)
            }
            let antiPure = max(anti[i].y - tail, 0.0)
            let expected = (pow(nuAS, 4) * nv) / (pow(nuS, 4) * (nv + 1.0))
            let measured = antiPure / aSi
            #expect(abs(measured / expected - 1.0) < 1e-3,
                    "ν_v=\(nuV)：anti/Stokes（去拖尾）=\(measured)，理论 \(expected)")
            #expect(abs(nv / (nv + 1.0) / exp(-hbar * omegaV / (kB * T)) - 1.0) < 1e-15,
                    "n_v/(n_v+1) 应等于玻尔兹曼因子 e^{−ħω/kT}")
            #expect(antiPure < aSi, "常温下 anti-Stokes 应弱于 Stokes")
        }
        // 金刚石 1332 cm⁻¹：室温布居 n_v ≈ 1.7×10⁻³（宝石学鉴别的一阶拉曼线）
        let diamond = try #require(db.ramanMode(id: "diamond_raman"))
        let omegaD = 2.0 * .pi * c_cm * diamond.nuVcm
        let nD = 1.0 / (exp(hbar * omegaD / (kB * T)) - 1.0)
        #expect(abs(nD - 1.685e-3) < 1e-5, "金刚石 1332 cm⁻¹ 室温布居 \(nD)")
    }

    @Test("物理律：anti-Stokes 随温度升高而增强（玻尔兹曼布居单调性）")
    func ramanTemperatureDependence() throws {
        let k = try constants()
        let hbar = try k.value("hbar"), kB = try k.value("kB"), c = try k.value("c")
        let c_cm = c * 100.0
        let nu0 = 1.0e7 / 785.0
        let mode = (nuV: 1332.0, fwhm: 6.0)
        let omegaV = 2.0 * .pi * c_cm * mode.nuV
        var previousNV = 0.0
        var previousRatio = 0.0
        for T in [50.0, 150.0, 300.0, 600.0] {
            let nv = 1.0 / (exp(hbar * omegaV / (kB * T)) - 1.0)
            let nuS = nu0 - mode.nuV, nuAS = nu0 + mode.nuV
            let ratio = (pow(nuAS, 4) * nv) / (pow(nuS, 4) * (nv + 1.0))
            #expect(nv > previousNV, "T=\(T)：布居 \(nv) 应高于上一档 \(previousNV)")
            #expect(ratio > previousRatio,
                    "T=\(T)：anti-Stokes/Stokes 比 \(ratio) 应高于上一档 \(previousRatio)")
            previousNV = nv
            previousRatio = ratio
        }
        // 高温极限 n_v ≫ 1 时比值 → (ν_AS/ν_S)⁴；低温极限 → 0
        // 注：1332 cm⁻¹ 模式在 600 K 时 n_v ≈ 4.3×10⁻²（仍 ≪ 1，但已显著非零）
        let nvLow = 1.0 / (exp(hbar * omegaV / (kB * 50.0)) - 1.0)
        let nvHigh = 1.0 / (exp(hbar * omegaV / (kB * 600.0)) - 1.0)
        #expect(nvLow < 1e-12 && nvHigh > 0.02, "低温几乎无 anti-Stokes，高温显著（n_v≈\(nvHigh)）")
    }

    // MARK: - fixtures 对拍

    @Test("fixtures：荧光/磷光衰减 400 点（混合容差 |Δ| ≤ max(1e-8, 1e-8·|y|)，尾部下溢为 0）")
    func decayMatchesFixture() async throws {
        let (_, tauH) = try hydrogenRates()
        let fx = try ModuleFixture.load("45_光谱学仪器的量子基础_荧光与寿命__v2022")
        let (tF, yF) = try fx.line(1, 0, index: 0)
        let (tP, yP) = try fx.line(1, 0, index: 1)
        #expect(tF.count == 400 && tP.count == 400)
        // t 网格对数均匀（1e-11 … 1e-2），fixture 9 位存储 → 相对核对
        #expect(abs(tF[0] / 1e-11 - 1.0) < 1e-8 && abs(tF[399] / 1e-2 - 1.0) < 1e-8)
        let cF = tF.map { exp(-$0 / tauH) }
        let cP = tP.map { exp(-$0 / 3.0e-3) }
        // 尾部（t ≫ τ）双精度下溢为 0：fixture 400 点中有 174 点为 0，相对误差在这些点失去意义，
        // 故统一用「绝对地板 1e-8 + 相对 1e-8」的混合判据（实测最大绝对偏差 2.0×10⁻⁹）。
        for i in 0..<400 {
            #expect(abs(cF[i] - yF[i]) <= max(1e-8, 1e-8 * abs(yF[i])),
                    "荧光第 \(i) 点：\(cF[i]) vs \(yF[i])")
            #expect(abs(cP[i] - yP[i]) <= max(1e-8, 1e-8 * abs(yP[i])),
                    "磷光第 \(i) 点：\(cP[i]) vs \(yP[i])")
        }
    }

    @Test("fixtures：拉曼谱 512 点（计算核求值 < 1e-6；60001 点重采样 < 2e-2，跳过峰芯 ±15 FWHM）")
    func ramanMatchesFixture() async throws {
        let k = try constants()
        let c = try k.value("c"), hbar = try k.value("hbar"), kB = try k.value("kB")
        let db = try MaterialDB.load()
        let modes = db.ramanModes.map { (nuV: $0.nuVcm, fwhm: $0.fwhmCm) }
        let fx = try ModuleFixture.load("45_光谱学仪器的量子基础_拉曼谱__v2022")
        let (x, y) = try fx.line(0, 0, index: 0)
        #expect(x.count == 512)
        let nu0 = 1.0e7 / 785.0
        #expect(abs(x[0] - (nu0 - 3000.0)) < 1e-3, "fixture 起点应为 ν₀−3000")
        let computed = SpectroscopyMath.ramanSpectrum(x, nu0: nu0, c_cm: c * 100.0,
                                                     hbar: hbar, kB: kB, T: 300.0,
                                                     modes: modes)
        // 容差 1e-6：fixture x 按 9 位有效数字存储（~1e-4 cm⁻¹），经 FWHM 6 cm⁻¹ 的
        // 洛伦兹斜率放大到 ~5×10⁻⁷（Python 独立复算实测最大相对偏差 5.3×10⁻⁷）。
        expectPointwiseClose(computed, y, tolerance: 1e-6, "拉曼谱（fixture 网格求值）")

        // 结构级重采样：fixture 512 点格距 11.7 cm⁻¹，宽于多数模式的线宽（6–14 cm⁻¹），
        // 峰芯无法插值（全量实测最大偏差 71%）；跳过峰芯 ±15 FWHM 后仅比对翼部与基线
        // （保留约 74% 样本），实测最大偏差 1.2%，取容差 2%。
        let module = RamanSpectrumModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        let chart = try #require(chartLineSeries(result, 0))
        let centers = modes.flatMap { [nu0 - $0.nuV, nu0 + $0.nuV] }
        let windows: [(Double, Double)] = modes.flatMap { m in
            [(nu0 - m.nuV, 15 * m.fwhm), (nu0 + m.nuV, 15 * m.fwhm)]
        }
        expectResampled(moduleX: chart.series[0].points.map(\.x),
                        moduleY: chart.series[0].points.map(\.y),
                        fx: x, fy: y, tolerance: 2e-2, "拉曼谱重采样",
                        skip: { xv in
                            windows.contains { abs(xv - $0.0) < $0.1 }
                        })
        #expect(centers.count == 6)
    }

    // MARK: - compute 结构与计时

    @Test("compute 结构：荧光寿命 3 图（双轴 300 点 + 衰减 400 点 + Jablonski 示意）+ 4 摘要")
    func fluorescenceStructureAndTiming() async throws {
        let module = FluorescenceLifetimeModule()
        let values = ParamValues.defaults(for: module.params)
        let result = try await module.compute(values, constants: try constants())
        #expect(result.charts.count == 3)
        let dual = try #require(chartDualAxis(result, 0))
        #expect(dual.primary[0].points.count == 300)
        #expect(dual.secondary[0].points.count == 300)
        let decay = try #require(chartLineSeries(result, 1))
        #expect(decay.series.count == 2)
        for s in decay.series { #expect(s.points.count == 400) }
        let jablonski = try #require(chartSchematic(result, 2))
        #expect(jablonski.callouts.count == 4, "吸收 / 荧光 / ISC / 磷光 四条标注")
        #expect(jablonski.referenceLines.count == 3, "S0 / S1 / T1 三条能级线")
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 3)
        try await expectComputeUnderBudget(module: module, values: values,
                                           constants: try constants(),
                                           budgetMillis: 16)
    }

    @Test("compute 结构：拉曼谱 1 图 60001 点 + 3 参考线 + 4 摘要（秒级档）")
    func ramanStructureAndTiming() async throws {
        let module = RamanSpectrumModule()
        let values = ParamValues.defaults(for: module.params)
        let result = try await module.compute(values, constants: try constants())
        #expect(result.charts.count == 1)
        let chart = try #require(chartLineSeries(result, 0))
        #expect(chart.series.count == 1)
        #expect(chart.series[0].points.count == 60001)
        #expect(chart.referenceLines.count == 3, "瑞利线 + 金刚石 Stokes/anti-Stokes")
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 4)
        // 玻尔兹曼一致性自检（模块摘要的 ✓/✗ 标记）
        let check = try #require(result.summary.first { $0.id == "Tcheck" })
        #expect(check.value == "✓", "n_v/(n_v+1) 与 e^{−ħω/kT} 应一致")
        try await expectComputeUnderBudget(module: module, values: values,
                                           constants: try constants(),
                                           budgetMillis: 2000)
    }
}
