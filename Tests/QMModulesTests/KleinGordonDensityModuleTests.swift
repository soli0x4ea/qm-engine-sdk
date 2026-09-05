import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W9-37 克莱因-戈登方程：概率密度振荡。
/// 定义式 ⇔ 解析式全网格恒等（脚本 assert < 1e-12）+ 双时刻剖面 fixture 逐点对拍
/// + 守恒律（空间平均不随 t）+ ρmin/ρmax/周期锚点。
@Suite("W9 克莱因-戈登方程：概率密度振荡")
struct KleinGordonDensityModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    // MARK: 计算核（物理律）

    @Test("物理律：定义式与解析式恒等（脚本口径全网格 max|Δ| < 1e-12）")
    func definitionMatchesAnalytic() {
        // 脚本默认参数 k1=1, k2=3, β=0.9（fixture 口径）
        let (k1, k2, beta) = (1.0, 3.0, 0.9)
        let lx = KleinGordonDensityMath.spatialPeriod(k1: k1, k2: k2)
        let lt = KleinGordonDensityMath.temporalPeriod(k1: k1, k2: k2)
        let xs = Num.linspace(0, lx, count: 50)
        let ts = Num.linspace(0, lt, count: 50)
        for t in ts {
            for x in xs {
                let byDef = KleinGordonDensityMath.density(x: x, t: t, k1: k1, k2: k2, beta: beta)
                let analytic = KleinGordonDensityMath.densityAnalytic(
                    x: x, t: t, k1: k1, k2: k2, beta: beta)
                #expect(abs(byDef - analytic) < 1e-12,
                        "x=\(x) t=\(t): \(byDef) vs \(analytic)")
            }
        }
    }

    @Test("物理律：空间平均 ∫ρdx 守恒（= ω₁ − ω₂β²，不随 t）")
    func spatialMeanConservation() {
        let (k1, k2, beta) = (1.0, 3.0, 0.9)
        let w1 = KleinGordonDensityMath.omega(k1)
        let w2 = KleinGordonDensityMath.omega(k2)
        let lx = KleinGordonDensityMath.spatialPeriod(k1: k1, k2: k2)
        let xs = Num.linspace(0, lx, count: 400)
        let expected = (w1 - w2 * beta * beta) * lx   // cos 项全周期积分为 0（∫ρdx，非平均）
        for frac in [0.0, 0.125, 0.25, 0.5, 0.75, 1.0] {
            let t = frac * KleinGordonDensityMath.temporalPeriod(k1: k1, k2: k2)
            let rho = xs.map {
                KleinGordonDensityMath.density(x: $0, t: t, k1: k1, k2: k2, beta: beta)
            }
            let integral = Num.trapezoid(xs, rho)
            #expect(abs(integral - expected) < 1e-9, "t/Lt = \(frac): ∫ρdx = \(integral)")
        }
    }

    @Test("脚本锚点：T = 3.5944；ρ_min = −2.7205；ρ_max = 0.4260；ρ>0 占比 0.240")
    func scriptAnchors() {
        let (k1, k2, beta) = (1.0, 3.0, 0.9)
        let lt = KleinGordonDensityMath.temporalPeriod(k1: k1, k2: k2)
        #expect(abs(lt - 3.5944) < 5e-4, "T = \(lt)")

        let lx = KleinGordonDensityMath.spatialPeriod(k1: k1, k2: k2)
        let xs = Num.linspace(0, lx, count: 400)
        let ts = Num.linspace(0, lt, count: 600)
        var rhoMin = Double.infinity, rhoMax = -Double.infinity
        var positive = 0
        for t in ts {
            for x in xs {
                let r = KleinGordonDensityMath.density(x: x, t: t, k1: k1, k2: k2, beta: beta)
                rhoMin = min(rhoMin, r)
                rhoMax = max(rhoMax, r)
                if r > 0 { positive += 1 }
            }
        }
        #expect(abs(rhoMin - (-2.7205)) < 5e-4, "ρ_min = \(rhoMin)")
        #expect(abs(rhoMax - 0.4260) < 5e-4, "ρ_max = \(rhoMax)")
        #expect(abs(Double(positive) / Double(400 * 600) - 0.240) < 5e-3)
        // 非正定：min < 0 < max
        #expect(rhoMin < 0 && rhoMax > 0)
    }

    // MARK: fixture 逐点对拍（双时刻剖面 400 点原生网格）

    @Test("fixture 对拍：t=0 与 t=T/2 空间剖面 400 点（rel < 1e-8）")
    func profileFixtures() throws {
        let fixture = try ModuleFixture.load("克莱因-戈登方程_概率密度振荡__v2022")
        let fx0 = try fixture.line(0, 0, label: "t=0")
        let fxHalf = try fixture.line(0, 0, label: "t=T/2")

        let (k1, k2, beta) = (1.0, 3.0, 0.9)
        let lt = KleinGordonDensityMath.temporalPeriod(k1: k1, k2: k2)
        let prof0 = fx0.x.map {
            KleinGordonDensityMath.density(x: $0, t: 0, k1: k1, k2: k2, beta: beta)
        }
        let profHalf = fxHalf.x.map {
            KleinGordonDensityMath.density(x: $0, t: lt / 2, k1: k1, k2: k2, beta: beta)
        }
        // fixture x 也存 9 位有效数字：零点附近 ρ 梯度 ~6.3 把 x 误差放大成
        // ~1e-8 的密度绝对误差（rel 可到 1e-5），用 abs + rel 组合判据
        func expectCombined(_ actual: [Double], _ expected: [Double], _ ctx: String) {
            #expect(actual.count == expected.count, "\(ctx): 点数不一致")
            for i in actual.indices {
                #expect(abs(actual[i] - expected[i]) <= 1e-7 + 1e-6 * abs(expected[i]),
                        "\(ctx) 第 \(i) 点：\(actual[i]) vs \(expected[i])")
            }
        }
        expectCombined(prof0, fx0.y, "t=0 剖面")
        expectCombined(profHalf, fxHalf.y, "t=T/2 剖面")
    }

    // MARK: 模块契约

    @Test("模块契约：等高线 400×600 + 剖面/时序 400 点 + 摘要 5 卡")
    func moduleContract() async throws {
        let module = KleinGordonDensityModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())

        #expect(result.charts.count == 3)
        guard case .contour(let contour) = result.charts[0] else {
            Issue.record("图 1 应为等高线"); return
        }
        #expect(contour.xGrid.count == 400 && contour.yGrid.count == 600)
        #expect(contour.values.count == 600 && contour.values[0].count == 400)

        guard case .lineSeries(let prof) = result.charts[1] else {
            Issue.record("图 2 应为曲线"); return
        }
        #expect(prof.series.count == 2)
        #expect(prof.series.allSatisfy { $0.points.count == 400 })
        guard case .lineSeries(let time) = result.charts[2] else {
            Issue.record("图 3 应为曲线"); return
        }
        #expect(time.series[0].points.count == 400)

        #expect(result.summary.count == 5)
        #expect(result.summary[2].value.contains("ρ_min"))
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("秒级档预算：compute < 2000 ms（400×600 解析网格）")
    func computeBudget() async throws {
        let module = KleinGordonDensityModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 2000, runs: 5, warmup: 2)
    }
}
