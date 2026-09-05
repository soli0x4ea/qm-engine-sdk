import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W7 笔记 42《晶体场与配位场理论》：d⁵ 高/低自旋交叉 / Tanabe-Sugano d² 能级图。
/// fixtures：42_晶体场与配位场理论_模型__v2022.json（400 点，与模块同网格）、
/// 42_晶体场与配位场理论_TanabeSugano__v2022.json（512 点，模块 600 点 → 重采样对拍）。
@Suite("W7 笔记42 晶体场与配位场理论")
struct CrystalFieldModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    // MARK: - 物理律：高/低自旋交叉

    @Test("物理律：自旋交叉点恰在 Δ_o = P——E_LS = 2(P − Δ_o) 的零点，两侧符号翻转")
    func spinCrossoverLaw() {
        let P = 20000.0
        #expect(CrystalFieldMath.highSpinEnergy(12345.0) == 0.0, "高自旋以 barycenter 为零（CFSE=0）")
        // 交叉点：E_LS(P) = E_HS = 0
        #expect(abs(CrystalFieldMath.lowSpinEnergy(delta: P, P: P)) < 1e-15)
        // Δ_o < P：低自旋高（高自旋稳定）；Δ_o > P：低自旋低（低自旋稳定）
        for delta in stride(from: 0.0, to: P, by: 1373.0) {
            #expect(CrystalFieldMath.lowSpinEnergy(delta: delta, P: P) > 0,
                    "Δ_o=\(delta) < P：高自旋占优")
        }
        for delta in stride(from: P + 1373.0, through: 40000.0, by: 1373.0) {
            #expect(CrystalFieldMath.lowSpinEnergy(delta: delta, P: P) < 0,
                    "Δ_o=\(delta) > P：低自旋占优")
        }
        // 斜率 dE_LS/dΔ_o = −2（CFSE = −20Dq = −2Δ_o 部分）
        let h = 1e-3
        let slope = (CrystalFieldMath.lowSpinEnergy(delta: 10000 + h, P: P)
                     - CrystalFieldMath.lowSpinEnergy(delta: 10000 - h, P: P)) / (2 * h)
        #expect(abs(slope + 2.0) < 1e-9, "dE_LS/dΔ_o = \(slope)")
        // Δ_o = 0 处 E_LS = 2P（两对成对能代价）
        #expect(abs(CrystalFieldMath.lowSpinEnergy(delta: 0, P: P) - 2 * P) < 1e-15)
        // 交叉点随 P 线性移动（判据的普适性）
        for pTest in [12000.0, 20000.0, 28500.0] {
            #expect(abs(CrystalFieldMath.lowSpinEnergy(delta: pTest, P: pTest)) < 1e-15)
        }
    }

    @Test("物理律：八面体 d 轨道分裂无迹（质心守恒）——1·12 + 3·2 + 3·(−6) = 0")
    func octahedralTraceFreeLaw() {
        // ^3F 在 O_h 场分裂（Dq 单位）：^3A2g=+12, ^3T2g=+2, ^3T1g(F)=−6
        // 简并度加权和（1:3:3）应为 0——配位场算符在 L 子空间内无迹
        for x in stride(from: 0.0, through: 40.0, by: 2.5) {
            let Dq = x / 10.0
            let trace = 1.0 * (12.0 * Dq) + 3.0 * (2.0 * Dq) + 3.0 * (-6.0 * Dq)
            #expect(abs(trace) < 1e-12, "Δ_o/B=\(x)：加权迹 \(trace)")
        }
        // 两支未混合能级之差恰为 (12 − 2) Dq = Δ_o，即图上横坐标本身：E(^3A2g) − E(^3T2g) = Δ_o/B
        for x in [0.0, 5.0, 12.5, 25.0, 40.0] {
            let lv = CrystalFieldMath.tanabeSuganoLevels(x: x)
            #expect(abs((lv.a2g - lv.t2g) - x) < 1e-12,
                    "Δ_o/B=\(x)：E(^3A2g)−E(^3T2g) = \(lv.a2g - lv.t2g)")
        }
    }

    @Test("物理律：^3T1g(F)–^3T1g(P) 的 2×2 久期方程——迹/行列式恒等 + 回避交叉间隙 ≥ 12√2")
    func secularEquationLaw() {
        // H = [[a, b], [b, c]]，b = −6√2（B 单位）；闭式 λ± = ((a+c) ± √((a−c)²+4b²))/2
        let b = -6.0 * sqrt(2.0)
        for (dTF, dTP) in [(-6.0, 15.0), (0.0, 15.0), (-24.0, 23.0), (-30.0, 10.0)] {
            let (lo, hi) = CrystalFieldMath.secularRoots(dT1F: dTF, dT1P: dTP)
            #expect(abs((lo + hi) - (dTF + dTP)) < 1e-12, "迹守恒：λ₋+λ₊ = a+c")
            #expect(abs((lo * hi) - (dTF * dTP - b * b)) < 1e-9, "行列式守恒：λ₋λ₊ = ac − b²")
            #expect(lo <= hi)
            // 回避交叉：两支间隙 ≥ 2|b| = 12√2，等号仅在 a = c 时成立
            #expect(hi - lo >= 2 * abs(b) - 1e-12, "间隙 \(hi - lo) < 12√2 = \(2 * abs(b))")
            #expect(abs((hi - lo) - sqrt((dTF - dTP) * (dTF - dTP) + 4 * b * b)) < 1e-12)
        }
        // 间隙极小值恰在 a = c（真交叉点，本模型在 Δ_o/B = −18.75，位于物理域外）
        let (l0, h0) = CrystalFieldMath.secularRoots(dT1F: 7.5, dT1P: 7.5)
        #expect(abs((h0 - l0) - 2 * abs(b)) < 1e-12, "a=c 时间隙取最小值 12√2")
    }

    @Test("物理律：Tanabe-Sugano 图 Δ_o/B = 25 处读数与公开图表一致（T2g≈23, A2g≈48, T1P≈36 B）")
    func tanabeSuganoReferenceReadings() {
        let ref = CrystalFieldMath.tanabeSuganoLevels(x: 25.0)
        // 公开 d² Tanabe-Sugano 图在 Δ_o/B=25 的读数（图表目视精度 ±1–3 B）
        #expect(abs(ref.t2g - 23.0) < 1.5, "^3T2g = \(ref.t2g) B（图表 ≈23）")
        #expect(abs(ref.a2g - 48.0) < 1.5, "^3A2g = \(ref.a2g) B（图表 ≈48）")
        #expect(abs(ref.t1p - 36.0) < 3.5, "^3T1g(P) = \(ref.t1p) B（图表 ≈36）")
        // 层序结构（d² 八面体，强场）：^3T1g(F) 基态 < ^3T2g < ^3T1g(P) < ^3A2g
        #expect(abs(ref.t1f) < 1e-15, "基态即能量零点")
        #expect(ref.t2g < ref.t1p && ref.t1p < ref.a2g,
                "层序 \(ref.t2g) < \(ref.t1p) < \(ref.a2g)")
        // 弱场极限 Δ_o → 0：^3T1g(P) 回到自由离子 ^3P 的 15 B（经耦合上移为 √(15²+288)）
        let zero = CrystalFieldMath.tanabeSuganoLevels(x: 0.0)
        #expect(abs(zero.t1p - sqrt(15.0 * 15.0 + 288.0)) < 1e-12,
                "Δ_o=0 时 ^3T1g(P) 相对基态 = √(15²+288) = \(sqrt(15.0 * 15.0 + 288.0))")
    }

    // MARK: - fixtures 对拍

    @Test("fixtures：高/低自旋双线 400 点（混合容差 |Δ| ≤ max(1e-3, 1e-8·|y|)，覆盖零点）")
    func crystalFieldModelMatchesFixture() async throws {
        let fx = try ModuleFixture.load("42_晶体场与配位场理论_模型__v2022")
        let P = 20000.0
        // line 0 = 高自旋（恒 0）；line 1 = 低自旋 2(P − Δ_o)；line 2 = 交叉竖线
        let (xHS, yHS) = try fx.line(0, 0, index: 0)
        let (xLS, yLS) = try fx.line(0, 0, index: 1)
        #expect(xHS.count == 400 && xLS.count == 400)
        let cross = try fx.line(0, 0, index: 2)
        #expect(cross.x == [20000.0, 20000.0], "交叉竖线在 Δ_o = P = 20000")
        #expect(yHS.allSatisfy { $0 == 0.0 }, "高自旋恒为 0")
        for i in xLS.indices {
            let c = CrystalFieldMath.lowSpinEnergy(delta: xLS[i], P: P)
            let f = yLS[i]
            // 低自旋线过零（Δ_o = P），故加 1e-3 绝对地板：fixture 9 位存储 × 斜率 −2 的放大
            #expect(abs(c - f) <= max(1e-3, 1e-8 * abs(f)),
                    "低自旋第 \(i) 点：\(c) vs \(f)")
        }
        // 模块曲线结构级对拍（同网格，重采样退化为同点插值）
        let module = CrystalFieldModelModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        let chart = try #require(chartLineSeries(result, 0))
        #expect(chart.series.count == 2)
        // 低自旋线在 Δ_o = P 处过零，插值残差（实测 ≤1e-4 cm⁻¹）在零点附近相对误差被放大，
        // 故跳过交叉点 ±1000 cm⁻¹，其余区间按相对容差 1e-6（实测 ≤3e-8）。
        expectResampled(moduleX: chart.series[1].points.map(\.x),
                        moduleY: chart.series[1].points.map(\.y),
                        fx: xLS, fy: yLS, tolerance: 1e-6, "低自旋 E_LS(Δ_o)",
                        skip: { abs($0 - 20000.0) < 1000.0 })
        let crossLine = try #require(chart.referenceLines.first)
        #expect(crossLine.value == 20000.0, "交叉参考线 Δ_o = P")
    }

    @Test("fixtures：Tanabe-Sugano 四支 512 点（计算核求值 < 1e-8；模块 600 点重采样 < 1e-5）")
    func tanabeSuganoMatchesFixture() async throws {
        let fx = try ModuleFixture.load("42_晶体场与配位场理论_TanabeSugano__v2022")
        let names = ["^3T1g(F)", "^3T2g(F)", "^3A2g(F)", "^3T1g(P)"]
        var fixtureCurves: [(x: [Double], y: [Double])] = []
        for i in 0..<4 {
            let (x, y) = try fx.line(0, 0, index: i)
            #expect(x.count == 512, "\(names[i]) 应为 512 点")
            fixtureCurves.append((x, y))
            let computed = x.map { CrystalFieldMath.tanabeSuganoLevels(x: $0) }
                .map { [$0.t1f, $0.t2g, $0.a2g, $0.t1p][i] }
            expectPointwiseClose(computed, y, tolerance: 1e-8, "TS \(names[i])")
        }

        let module = CrystalFieldTanabeSuganoModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        let chart = try #require(chartLineSeries(result, 0))
        #expect(chart.series.count == 4)
        for (i, series) in chart.series.enumerated() {
            #expect(series.points.count == 600)
            expectResampled(moduleX: series.points.map(\.x), moduleY: series.points.map(\.y),
                            fx: fixtureCurves[i].x, fy: fixtureCurves[i].y,
                            tolerance: 1e-5, "TS \(names[i]) 重采样")
        }
    }

    // MARK: - compute 结构与计时

    @Test("compute 结构：高/低自旋 1 图双系列 400 点 + 交叉参考线 + 3 摘要")
    func crystalFieldModelStructureAndTiming() async throws {
        let module = CrystalFieldModelModule()
        let values = ParamValues.defaults(for: module.params)
        let result = try await module.compute(values, constants: try constants())
        #expect(result.charts.count == 1)
        let chart = try #require(chartLineSeries(result, 0))
        #expect(chart.series.count == 2)
        for s in chart.series { #expect(s.points.count == 400) }
        #expect(chart.referenceLines.count == 1)
        #expect(result.summary.count == 3)
        #expect(result.theory?.formulas.count == 3)
        // 换 P 值后交叉点随之移动（参数联动）
        var shifted = values
        shifted.sliders["P"] = 12000
        let r2 = try await module.compute(shifted, constants: try constants())
        let chart2 = try #require(chartLineSeries(r2, 0))
        let cross2 = try #require(chart2.referenceLines.first)
        #expect(cross2.value == 12000, "P 改为 12000 后交叉点同步移动")
        try await expectComputeUnderBudget(module: module, values: values,
                                           constants: try constants(),
                                           budgetMillis: 16)
    }

    @Test("compute 结构：Tanabe-Sugano 1 图四系列 600 点 + 对照参考线 + 4 摘要（秒级档）")
    func tanabeSuganoStructureAndTiming() async throws {
        let module = CrystalFieldTanabeSuganoModule()
        let values = ParamValues.defaults(for: module.params)
        let result = try await module.compute(values, constants: try constants())
        #expect(result.charts.count == 1)
        let chart = try #require(chartLineSeries(result, 0))
        #expect(chart.series.count == 4)
        for s in chart.series { #expect(s.points.count == 600) }
        #expect(chart.referenceLines.count == 1)
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 4)
        try await expectComputeUnderBudget(module: module, values: values,
                                           constants: try constants(),
                                           budgetMillis: 2000)
    }
}
