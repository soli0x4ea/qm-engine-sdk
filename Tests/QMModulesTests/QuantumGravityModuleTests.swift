import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W7 笔记 41《量子引力与全息原理》：霍金温度 / Bekenstein-Hawking 黑洞熵。
/// fixtures：量子引力与全息原理_霍金温度__v2022.json（log-log 400 点，M 从 m_P 到 1e10 M☉）、
/// 量子引力与全息原理_黑洞熵__v2022.json（400 点双曲线：面积律 + 普朗克单位）。
///
/// 对拍口径（与 W5/W6 一致）：fixtures 的 x/y 均按 9 位有效数字存储，引入 ~1e-8 相对差，
/// 故判据级对照一律「把计算核在 fixture 自身 x 网格上求值」（容差 1e-8， houses 惯例；
/// 计划 §6.1 一档 1e-12 是相对 Python 全精度的口径，受 fixture 存储精度封顶实际取 1e-8），
/// 结构级对照再用 `expectResampled` 比对模块 compute 产出的曲线。
@Suite("W7 笔记41 量子引力与全息原理")
struct QuantumGravityModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    private struct Phys {
        let hbar: Double, c: Double, G: Double, kB: Double
        let mP: Double, lP: Double
    }

    private func phys() throws -> Phys {
        let k = try constants()
        let hbar = try k.value("hbar"), c = try k.value("c")
        let G = try k.value("G"), kB = try k.value("kB")
        return Phys(hbar: hbar, c: c, G: G, kB: kB,
                    mP: QuantumGravityMath.planckMass(hbar: hbar, c: c, G: G),
                    lP: QuantumGravityMath.planckLength(hbar: hbar, c: c, G: G))
    }

    // MARK: - 物理律：霍金温度

    @Test("物理律：T_H ∝ M⁻¹——log-log 斜率恰为 −1，且半质量翻倍温度减半")
    func hawkingInverseMassLaw() async throws {
        let p = try phys()
        let module = HawkingTemperatureModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        let series = try #require(chartLineSeries(result, 0)?.series.first)
        #expect(series.points.count == 400)
        let first = try #require(series.points.first)
        let last = try #require(series.points.last)
        // 端点连线斜率 = −1（T_H ∝ M⁻¹，log-log 下为直线）
        let slope = (log10(last.y) - log10(first.y)) / (log10(last.x) - log10(first.x))
        #expect(abs(slope + 1.0) < 1e-12, "实测斜率 \(slope)")
        // 逐点斜率恒为 −1（闭合形式，非仅端点）
        for i in 1..<series.points.count {
            let a = series.points[i - 1], b = series.points[i]
            let s = (log10(b.y) - log10(a.y)) / (log10(b.x) - log10(a.x))
            #expect(abs(s + 1.0) < 1e-12, "第 \(i) 段斜率 \(s)")
        }
        // 质量加倍 → 温度减半（严格 2 倍反比）
        for M in [1.0e20, 1.98847e30, 6.5e9 * 1.98847e30] {
            let t1 = QuantumGravityMath.hawkingTemperature(M: M, hbar: p.hbar, c: p.c, G: p.G, kB: p.kB)
            let t2 = QuantumGravityMath.hawkingTemperature(M: 2 * M, hbar: p.hbar, c: p.c, G: p.G, kB: p.kB)
            #expect(abs(t1 / t2 - 2.0) < 1e-15, "M=\(M)：T(M)/T(2M) 应为 2")
        }
    }

    @Test("物理律：绝对量级核对——m_P≈2.176×10⁻⁸ kg、T_H(1 M☉)≈6.17×10⁻⁸ K（远低于 CMB）")
    func hawkingAbsoluteScale() async throws {
        let p = try phys()
        // 普朗克质量 2.176434×10⁻⁸ kg（CODATA 2018/2022 同值）
        #expect(abs(p.mP / 2.176434e-8 - 1.0) < 1e-6, "m_P = \(p.mP) kg")
        let tSun = QuantumGravityMath.hawkingTemperature(
            M: QuantumGravityMath.mSun, hbar: p.hbar, c: p.c, G: p.G, kB: p.kB)
        #expect(abs(tSun / 6.171e-8 - 1.0) < 2e-3, "T_H(1 M☉) = \(tSun) K")
        #expect(tSun < 2.725, "恒星级黑洞比 CMB 冷（净吸积而非蒸发）")
        // M87*（6.5×10⁹ M☉）温度 = 太阳质量黑洞 / 6.5×10⁹
        let tM87 = QuantumGravityMath.hawkingTemperature(
            M: QuantumGravityMath.m87, hbar: p.hbar, c: p.c, G: p.G, kB: p.kB)
        #expect(abs(tSun / tM87 / 6.5e9 - 1.0) < 1e-15)
        // 史瓦西半径线性于 M：R_s(2M) = 2 R_s(M)
        let r1 = QuantumGravityMath.schwarzschildRadius(M: 1.98847e30, c: p.c, G: p.G)
        let r2 = QuantumGravityMath.schwarzschildRadius(M: 2 * 1.98847e30, c: p.c, G: p.G)
        #expect(abs(r2 / r1 - 2.0) < 1e-15)
        #expect(abs(r1 / 2953.25 - 1.0) < 1e-3, "R_s(1 M☉) ≈ 2.953 km")
    }

    // MARK: - 物理律：黑洞熵

    @Test("物理律：面积律 S = A/(4 l_P²) 与普朗克单位 4π(M/m_P)² 完全等价（面积律 vs 体积律）")
    func entropyAreaLaw() async throws {
        let p = try phys()
        for M in [1.0e12, 1.98847e30, 1.0e36, QuantumGravityMath.m87] {
            let sArea = QuantumGravityMath.bekensteinHawkingS(M: M, hbar: p.hbar, c: p.c, G: p.G)
            let sPlanck = 4 * .pi * pow(M / p.mP, 2)
            #expect(abs(sArea / sPlanck - 1.0) < 1e-12,
                    "M=\(M)：面积律 \(sArea) vs 普朗克单位 \(sPlanck)")
            // S·(4 l_P²) 恰为视界面积 A = 4π R_s²
            let rs = QuantumGravityMath.schwarzschildRadius(M: M, c: p.c, G: p.G)
            let area = 4 * .pi * rs * rs
            #expect(abs(sArea * 4 * p.lP * p.lP / area - 1.0) < 1e-12)
        }
        // 普朗克质量黑洞的熵 = 4π（全纯数值，无单位残留）
        #expect(abs(QuantumGravityMath.bekensteinHawkingS(
            M: p.mP, hbar: p.hbar, c: p.c, G: p.G) / (4 * .pi) - 1.0) < 1e-12)
    }

    @Test("物理律：S ∝ M²（log-log 斜率 +2）——熵随质量平方暴涨，正比面积而非体积")
    func entropyMassSquaredLaw() async throws {
        let p = try phys()
        let module = BlackHoleEntropyModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        let chart = try #require(chartLineSeries(result, 0))
        #expect(chart.series.count == 2)
        let main = chart.series[0]
        #expect(main.points.count == 400)
        for i in 1..<main.points.count {
            let a = main.points[i - 1], b = main.points[i]
            let s = (log10(b.y) - log10(a.y)) / (log10(b.x) - log10(a.x))
            #expect(abs(s - 2.0) < 1e-12, "第 \(i) 段斜率 \(s)")
        }
        // 第二支（普朗克单位，M 从 m_P 起）在重叠区间与第一支完全重合：
        // 两支 x 网格不同（起点 m_P vs 1 M☉），故按闭式 4π(M/m_P)² 逐点核对第一支。
        let alt = chart.series[1]
        #expect(alt.points.count == 400)
        let altFirst = try #require(alt.points.first)
        #expect(altFirst.x < main.points[0].x, "普朗克单位支起点应更小（m_P < 1 M☉）")
        for pt in main.points {
            let closed = 4 * .pi * pow(pt.x / p.mP, 2)
            #expect(abs(pt.y / closed - 1.0) < 1e-9, "M=\(pt.x)：\(pt.y) vs \(closed)")
        }
        // 质量 ×10 → 熵 ×100
        let s1 = QuantumGravityMath.bekensteinHawkingS(M: 1.98847e30, hbar: p.hbar, c: p.c, G: p.G)
        let s2 = QuantumGravityMath.bekensteinHawkingS(M: 1.98847e31, hbar: p.hbar, c: p.c, G: p.G)
        #expect(abs(s2 / s1 - 100.0) < 1e-12)
    }

    // MARK: - fixtures 对拍

    @Test("fixtures：霍金温度 log-log 400 点（计算核在 fixture 网格求值 < 1e-12；重采样 < 1e-8）")
    func hawkingMatchesFixture() async throws {
        let p = try phys()
        let fx = try ModuleFixture.load("量子引力与全息原理_霍金温度__v2022")
        let (x, y) = try fx.line(0, 0, index: 0)
        #expect(x.count == 400)
        let computed = x.map {
            QuantumGravityMath.hawkingTemperature(M: $0, hbar: p.hbar, c: p.c, G: p.G, kB: p.kB)
        }
        expectPointwiseClose(computed, y, tolerance: 1e-8, "T_H(M)")

        let module = HawkingTemperatureModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        let series = try #require(chartLineSeries(result, 0)?.series.first)
        expectResampled(moduleX: series.points.map(\.x), moduleY: series.points.map(\.y),
                        fx: x, fy: y, tolerance: 1e-8, "霍金温度 log-log", logSpace: true)
    }

    @Test("fixtures：黑洞熵双曲线 400 点（面积律 + 普朗克单位，各自 < 1e-12 / 1e-8）")
    func entropyMatchesFixture() async throws {
        let p = try phys()
        let fx = try ModuleFixture.load("量子引力与全息原理_黑洞熵__v2022")
        // line 0 = 面积律（M ≥ 1 M☉）；line 1 = 普朗克单位（M ≥ m_P）
        let (xArea, yArea) = try fx.line(0, 0, index: 0)
        let (xPlanck, yPlanck) = try fx.line(0, 0, index: 1)
        #expect(xArea.count == 400 && xPlanck.count == 400)
        let cArea = xArea.map { QuantumGravityMath.bekensteinHawkingS(M: $0, hbar: p.hbar, c: p.c, G: p.G) }
        expectPointwiseClose(cArea, yArea, tolerance: 1e-8, "S_BH(M) 面积律")
        let cPlanck = xPlanck.map { 4 * .pi * pow($0 / p.mP, 2) }
        expectPointwiseClose(cPlanck, yPlanck, tolerance: 1e-8, "S_BH(M) 普朗克单位")

        let module = BlackHoleEntropyModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        let chart = try #require(chartLineSeries(result, 0))
        expectResampled(moduleX: chart.series[0].points.map(\.x),
                        moduleY: chart.series[0].points.map(\.y),
                        fx: xArea, fy: yArea, tolerance: 1e-8, "黑洞熵·面积律", logSpace: true)
        expectResampled(moduleX: chart.series[1].points.map(\.x),
                        moduleY: chart.series[1].points.map(\.y),
                        fx: xPlanck, fy: yPlanck, tolerance: 1e-8, "黑洞熵·普朗克单位", logSpace: true)
    }

    // MARK: - compute 结构与计时

    @Test("compute 结构：霍金温度 1 图 400 点 + 4 参考线 + 4 摘要 + 理论卡")
    func hawkingStructureAndTiming() async throws {
        let module = HawkingTemperatureModule()
        let values = ParamValues.defaults(for: module.params)
        let result = try await module.compute(values, constants: try constants())
        #expect(result.charts.count == 1)
        let chart = try #require(chartLineSeries(result, 0))
        #expect(chart.series.count == 1)
        #expect(chart.series[0].points.count == 400)
        #expect(chart.referenceLines.count == 4, "m_P / 1 M☉ / 10 M☉ / M87* 四条参考线")
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 3)
        try await expectComputeUnderBudget(module: module, values: values,
                                           constants: try constants(),
                                           budgetMillis: 16)
    }

    @Test("compute 结构：黑洞熵 1 图双系列各 400 点 + 2 参考线 + 4 摘要 + 理论卡")
    func entropyStructureAndTiming() async throws {
        let module = BlackHoleEntropyModule()
        let values = ParamValues.defaults(for: module.params)
        let result = try await module.compute(values, constants: try constants())
        #expect(result.charts.count == 1)
        let chart = try #require(chartLineSeries(result, 0))
        #expect(chart.series.count == 2)
        for s in chart.series { #expect(s.points.count == 400) }
        #expect(chart.referenceLines.count == 2)
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 3)
        try await expectComputeUnderBudget(module: module, values: values,
                                           constants: try constants(),
                                           budgetMillis: 16)
    }
}
