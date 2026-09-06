import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W7 笔记 43《过渡金属离子致色的量子机制》：振动耦合饱和强度 / λ = 10⁷/Δ 致色模型。
/// fixtures：43_过渡金属离子致色的量子机制_振动耦合__v2022.json（512 点，模块 600 点）、
/// 43_过渡金属离子致色的量子机制_模型__v2022.json（400 点同网格 + 3 个矿物散点）。
@Suite("W7 笔记43 过渡金属离子致色的量子机制")
struct ColorMechanismModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    // MARK: - 物理律：振动耦合饱和

    @Test("物理律：I(ξ) = ξ²/(ξ²+ξ₀²)——半强度点、饱和上限、倒易对称 I(ξ)+I(ξ₀²/ξ)=1")
    func vibronicSaturationLaw() {
        let xi0 = 0.5
        #expect(ColorMechanismMath.relativeIntensity(xi: 0, xi0: xi0) == 0)
        #expect(abs(ColorMechanismMath.relativeIntensity(xi: xi0, xi0: xi0) - 0.5) < 1e-15,
                "ξ = ξ₀ 处恰为半强度")
        // 单调性
        for xi in stride(from: 0.0, to: 3.0, by: 0.05) {
            let a = ColorMechanismMath.relativeIntensity(xi: xi, xi0: xi0)
            let b = ColorMechanismMath.relativeIntensity(xi: xi + 0.05, xi0: xi0)
            #expect(b > a, "I 应严格单调递增（ξ=\(xi)）")
        }
        // 倒易对称：I(ξ) + I(ξ₀²/ξ) = 1（饱和函数的配对性质）
        for xi in stride(from: 0.05, through: 3.0, by: 0.05) {
            let s = ColorMechanismMath.relativeIntensity(xi: xi, xi0: xi0)
                + ColorMechanismMath.relativeIntensity(xi: xi0 * xi0 / xi, xi0: xi0)
            #expect(abs(s - 1.0) < 1e-15, "ξ=\(xi)：配对和 \(s)")
        }
        // 饱和上限 1（ξ → ∞）：ξ = 100 ξ₀ 时已到 0.9999
        let big = ColorMechanismMath.relativeIntensity(xi: 100 * xi0, xi0: xi0)
        #expect(abs(big - 10000.0 / 10001.0) < 1e-12)
        #expect(big < 1.0, "相对强度恒 < 1（归一化）")
        // 小 ξ 极限：I = (ξ/ξ₀)²/(1 + (ξ/ξ₀)²) —— 对任意 ξ 都是精确恒等式
        for xi in [1e-3, 1e-2, 0.1, 0.5, 1.0, 3.0] {
            let r = xi * xi / (xi0 * xi0)
            let ratio = ColorMechanismMath.relativeIntensity(xi: xi, xi0: xi0) / r
            #expect(abs(ratio - 1.0 / (1.0 + r)) < 1e-15, "ξ=\(xi)：I/(ξ/ξ₀)² = \(ratio)")
        }
    }

    @Test("物理律：梯形积分 ∫₀³ I(ξ)dξ = 3 − ξ₀·arctan(3/ξ₀)（解析值 2.29718，ξ₀=0.5）")
    func vibronicAreaLaw() async throws {
        let xi0 = 0.5
        let xi = Num.linspace(0.0, 3.0, count: 600)
        let I = xi.map { ColorMechanismMath.relativeIntensity(xi: $0, xi0: xi0) }
        let area = Num.trapezoid(xi, I)
        let analytic = 3.0 - xi0 * atan(3.0 / xi0)
        #expect(abs(area - analytic) < 1e-6, "数值积分 \(area) vs 解析 \(analytic)")
        #expect(abs(analytic - 2.2971762) < 1e-6, "解析值 ≈ 2.29718")
        // 模块 compute 摘要里的面积卡应与解析值一致（%.3f 格式下为 2.297）
        let module = ColorVibrationalCouplingModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        let areaCard = try #require(result.summary.first { $0.id == "area" })
        #expect(areaCard.value == String(format: "%.3f", analytic),
                "摘要面积 \(areaCard.value) 应等于解析值 \(analytic)")
    }

    // MARK: - 物理律：λ = 10⁷/Δ 反比律

    @Test("物理律：λ·Δ = 10⁷ nm·cm⁻¹（log-log 斜率 −1），MaterialDB 三行矿物均落在曲线上")
    func inverseWavelengthLaw() throws {
        // 曲线上的反比律
        for delta in stride(from: 5000.0, through: 25000.0, by: 1000.0) {
            let lam = ColorMechanismMath.deltaToLambdaNm(deltaCm: delta)
            #expect(abs(lam * delta - 1.0e7) < 1e-6, "Δ=\(delta)：λ·Δ = \(lam * delta)")
        }
        let a = ColorMechanismMath.deltaToLambdaNm(deltaCm: 8000)
        let b = ColorMechanismMath.deltaToLambdaNm(deltaCm: 16000)
        #expect(abs(a / b - 2.0) < 1e-15, "Δ 减半 → λ 加倍")

        // MaterialDB 三行矿物必须落在同一条 λ = 10⁷/Δ 曲线上
        let db = try MaterialDB.load()
        #expect(db.crystalFieldMinerals.count == 3)
        var previousDelta = Double.infinity
        var previousLambda = 0.0
        for m in db.crystalFieldMinerals {
            let lam = ColorMechanismMath.deltaToLambdaNm(deltaCm: m.deltaCm)
            #expect(abs(lam * m.deltaCm - 1.0e7) < 1e-6, "\(m.id) 不在 λ=10⁷/Δ 曲线上")
            // Δ_o 递减 ↔ λ 递增（同一 Cr³⁺ 在刚玉中偏红、在绿柱石中偏绿）
            #expect(m.deltaCm < previousDelta, "\(m.id) 的 Δ_o 未按降序登记")
            #expect(lam > previousLambda, "\(m.id) 的 λ 应按升序")
            previousDelta = m.deltaCm
            previousLambda = lam
        }
        // 红宝石/祖母绿在可见光区（400–780 nm），橄榄石主带在近红外
        let ruby = try #require(db.mineral(id: "ruby"))
        let emerald = try #require(db.mineral(id: "emerald"))
        let peridot = try #require(db.mineral(id: "peridot"))
        #expect(ColorMechanismMath.deltaToLambdaNm(deltaCm: ruby.deltaCm) < 780)
        #expect(ColorMechanismMath.deltaToLambdaNm(deltaCm: emerald.deltaCm) < 780)
        #expect(ColorMechanismMath.deltaToLambdaNm(deltaCm: peridot.deltaCm) > 780,
                "橄榄石 Fe²⁺ 主带在近红外（~1050 nm）")
    }

    @Test("物理律：O_h 八面体分裂——质心守恒 2E(eg)+3E(t2g)=0，Δ_o = 10Dq，能级比 −3:2")
    func octahedralSplittingLaw() {
        for dq in [500.0, 1830.0, 3000.0] {
            let sp = ColorMechanismMath.octahedralSplitting(dqCm: dq)
            #expect(abs(sp.barycenter) < 1e-12, "质心守恒：2E(eg)+3E(t2g) = \(sp.barycenter)")
            #expect(abs(sp.delta - 10 * dq) < 1e-12, "Δ_o = 10Dq")
            #expect(abs(sp.eEg / sp.eT2g + 1.5) < 1e-15, "E(eg)/E(t2g) = 6/(−4) = −3/2")
            #expect(abs(sp.eEg - 6 * dq) < 1e-12 && abs(sp.eT2g + 4 * dq) < 1e-12)
        }
        // 可见光窗 Δ ∈ [12820, 25000] cm⁻¹ ↔ λ ∈ [400, 780] nm
        #expect(abs(ColorMechanismMath.deltaToLambdaNm(deltaCm: 12820) - 780.03) < 0.05)
        #expect(abs(ColorMechanismMath.deltaToLambdaNm(deltaCm: 25000) - 400.0) < 1e-9)
    }

    // MARK: - fixtures 对拍

    @Test("fixtures：振动耦合饱和曲线 512 点（计算核求值 < 1e-8；重采样 < 1e-3，跳过 ξ<0.1）")
    func vibronicMatchesFixture() async throws {
        let fx = try ModuleFixture.load("43_过渡金属离子致色的量子机制_振动耦合__v2022")
        let (x, y) = try fx.line(0, 0, index: 0)
        #expect(x.count == 512)
        let xi0 = 0.5
        let computed = x.map { ColorMechanismMath.relativeIntensity(xi: $0, xi0: xi0) }
        expectPointwiseClose(computed, y, tolerance: 1e-8, "I(ξ)")

        let module = ColorVibrationalCouplingModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        let chart = try #require(chartLineSeries(result, 0))
        #expect(chart.series.count == 1)
        // W13h #43：显示窗随 ξ₀ 收缩（默认 ξ₀=0.5 → x ≤ 2.0 → ~400 点）
        #expect(chart.series[0].points.count >= 350 && chart.series[0].points.count <= 450,
                "窗内点数 \(chart.series[0].points.count)")
        // ξ → 0 段 I ∝ ξ²，512 点 fixture 的线性插值在二次起始段相对误差被放大到 10%，
        // 故跳过 ξ < 0.1（I < 4%）；其余区间实测插值误差 ≤ 1.3×10⁻³，取容差 5×10⁻³。
        expectResampled(moduleX: chart.series[0].points.map(\.x),
                        moduleY: chart.series[0].points.map(\.y),
                        fx: x, fy: y, tolerance: 5e-3, "I(ξ) 重采样",
                        skip: { $0 < 0.1 })
    }

    @Test("fixtures：λ = 10⁷/Δ 曲线 400 点 + 三矿物散点（MaterialDB 与 fixture 散点一致）")
    func colorModelMatchesFixture() async throws {
        let fx = try ModuleFixture.load("43_过渡金属离子致色的量子机制_模型__v2022")
        let (x, y) = try fx.line(0, 0, index: 0)
        #expect(x.count == 400)
        let computed = x.map { ColorMechanismMath.deltaToLambdaNm(deltaCm: $0) }
        expectPointwiseClose(computed, y, tolerance: 1e-8, "λ(Δ)")

        let scatters = try fx.scatters(0, 0)
        #expect(scatters.count == 3, "fixture 应有 3 个矿物散点")

        let module = ColorMechanismModelModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        let curve = try #require(chartLineSeries(result, 0))
        #expect(curve.series[0].points.count == 400)
        expectResampled(moduleX: curve.series[0].points.map(\.x),
                        moduleY: curve.series[0].points.map(\.y),
                        fx: x, fy: y, tolerance: 1e-4, "λ(Δ) 重采样")

        // 散点：模块（MaterialDB）↔ fixture（Python 管线录入的同一批数据）
        let scatter = try #require(chartScatter(result, 1))
        #expect(scatter.series.count == 3)
        for (i, series) in scatter.series.enumerated() {
            let pt = try #require(series.points.first)
            #expect(series.points.count == 1)
            #expect(abs(pt.x / scatters[i].x[0] - 1.0) < 1e-8,
                    "散点 \(i) 的 Δ：\(pt.x) vs \(scatters[i].x[0])")
            #expect(abs(pt.y / scatters[i].y[0] - 1.0) < 1e-7,
                    "散点 \(i) 的 λ：\(pt.y) vs \(scatters[i].y[0])")
        }
    }

    // MARK: - compute 结构与计时

    @Test("compute 结构：振动耦合 1 图（显示窗随 ξ₀ 收缩，W13h）+ 双参考线 + 4 摘要")
    func vibronicStructureAndTiming() async throws {
        let module = ColorVibrationalCouplingModule()
        let values = ParamValues.defaults(for: module.params)
        let result = try await module.compute(values, constants: try constants())
        #expect(result.charts.count == 1)
        let chart = try #require(chartLineSeries(result, 0))
        #expect(chart.series.count == 1)
        // W13h #43：ξ₀=0.5 → 显示窗 x ≤ min(3, max(0.6, 4ξ₀)) = 2.0 → 600·(2/3) ≈ 400 点
        let pts = chart.series[0].points
        #expect(pts.count >= 350 && pts.count <= 450, "窗内点数 \(pts.count)")
        #expect(pts.allSatisfy { $0.x <= 2.0 + 1e-9 }, "默认窗 x ≤ 2.0")
        // ξ₀ 竖线 + 端点标记在窗内可见
        #expect(chart.referenceLines.count == 2)
        #expect(chart.pointMarkers.count == 1)
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 3)
        try await expectComputeUnderBudget(module: module, values: values,
                                           constants: try constants(),
                                           budgetMillis: 16)
    }

    @Test("compute 结构：致色模型 2 图（400 点曲线 + 3 散点）+ 可见光双参考线 + 10Dq 线 + 6 摘要")
    func colorModelStructureAndTiming() async throws {
        let module = ColorMechanismModelModule()
        let values = ParamValues.defaults(for: module.params)
        let result = try await module.compute(values, constants: try constants())
        #expect(result.charts.count == 2)
        let curve = try #require(chartLineSeries(result, 0))
        #expect(curve.series[0].points.count == 400)
        #expect(curve.referenceLines.count == 3, "可见光上下界 + 10Dq（W13h）")
        #expect(curve.pointMarkers.count == 1, "W13h：10Dq → λ 标记点")
        #expect(abs(curve.pointMarkers[0].x - 10 * values.sliders["Dq"]!) < 1e-6)
        let scatter = try #require(chartScatter(result, 1))
        #expect(scatter.series.count == 3)
        #expect(result.summary.count == 6, "3 矿物 + 反比律 + 分裂 + 质心")
        #expect(result.theory?.formulas.count == 3)
        try await expectComputeUnderBudget(module: module, values: values,
                                           constants: try constants(),
                                           budgetMillis: 16)
    }
}
