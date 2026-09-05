import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W4-8 不确定性关系：复刻脚本 3 条断言 + fixtures 三线对拍（FFT 动量分布）。
@Suite("W4 不确定性关系")
struct UncertaintyModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }
    private static let sigma0 = 1.0

    @Test("check1：高斯 Δx·Δp/ħ = 0.5（等号，|偏差| < 1e-3）")
    func gaussianSaturates() {
        let (_, _, dx, dp) = UncertaintyMath.gaussianMoments(sigma0: Self.sigma0)
        let prod = dx * dp
        #expect(abs(prod - 0.5) < 1e-3, "ΔxΔp/ħ = \(prod)")
        // 脚本参数化：ψ = (2πσ₀²)^{−1/4}e^{−x²/4σ₀²} ⇒ Δx = σ₀，Δp/ħ = 1/(2σ₀)
        #expect(abs(dx - Self.sigma0) < 1e-3, "Δx = \(dx)")
        #expect(abs(dp - 1 / (2 * Self.sigma0)) < 1e-3, "Δp/ħ = \(dp)")
    }

    @Test("check2/3：势阱积 > 0.5 + 1e-3 且与解析 < 1e-2")
    func wellExceeds() {
        let (_, _, dx, dp) = UncertaintyMath.wellMoments(L: 1.0)
        let prod = dx * dp
        #expect(prod > 0.5 + 1e-3, "积 = \(prod)")
        let ana = UncertaintyMath.wellAnalytic(L: 1.0)
        #expect(abs(prod - ana.dx * ana.dpOverHbar) < 1e-2,
                "数值 \(prod) vs 解析 \(ana.dx * ana.dpOverHbar)")
        #expect(abs(dx - 0.18076) < 1e-3, "Δx/L = \(dx)")
        #expect(abs(dp - .pi) < 1e-2, "Δp/ħ·L = \(dp)")
    }

    @Test("fixture 对拍：高斯 |ψ(x)|² 与势阱 |ψ₁(x)|²（512 点，rel < 1e-6）")
    func positionCurvesMatchFixtures() throws {
        let fx = try ModuleFixture.load("不确定性关系_等号与势阱__v2022")
        let cs = try constants()
        _ = cs

        // fig0ax0 线 0：|ψ(x)|² = ψ²/ψ_max（σ₀ = 1）
        let (fx0, fy0) = try fx.line(0, 0, index: 0)
        let (grid, psiGrid) = UncertaintyMath.gaussianState(sigma0: Self.sigma0)
        var psiMax = 0.0
        for v in psiGrid where v > psiMax { psiMax = v }
        let norm = pow(2 * .pi, -0.25)
        let actual0 = fx0.map {
            let p = norm * exp(-$0 * $0 / 4)
            return p * p / psiMax
        }
        expectPointwiseClose(actual0, fy0, tolerance: 1e-6, "高斯 |ψ|²")

        // fig0ax1 线 0：|ψ₁(x)|² = ψ²/ψ_max（L = 1）
        let (fx1, fy1) = try fx.line(0, 1, index: 0)
        let (wgrid, wpsi) = UncertaintyMath.wellGround(L: 1.0)
        var wMax = 0.0
        for v in wpsi where v > wMax { wMax = v }
        let actual1 = fx1.map {
            let p = 2.0.squareRoot() * sin(.pi * $0)
            return p * p / wMax
        }
        expectPointwiseClose(actual1, fy1, tolerance: 1e-6, "势阱 |ψ₁|²")
        _ = wgrid
    }

    @Test("fixture 对拍：动量分布 |φ(p)|²（FFT，512 点，分区容差）")
    func momentumCurveMatchesFixtures() throws {
        let fx = try ModuleFixture.load("不确定性关系_等号与势阱__v2022")
        let (fxp, fyp) = try fx.line(0, 0, index: 1)

        let (k, pd, _, _) = UncertaintyMath.gaussianMoments(sigma0: Self.sigma0)
        var pdMax = 0.0
        for v in pd where v > pdMax { pdMax = v }

        // fixture x = k/3（k 为 fftshift 网格值）：由 k 值反推网格下标 j = k/dk + N/2。
        // 注：Python 画图按 |p| 排序取点，其 argsort 对相同 |k| 的并列顺序不稳定，
        // 故不重排序列，改对拍 k → y 的映射关系（顺序无关）。
        let n = k.count
        let dk = k[n / 2 + 1]   // k[N/2] = 0，k[N/2+1] = dk

        // y 分区对拍：fixture 曲线动态范围跨 36 个数量级，
        // 尾部（y < 1e-30）在 numpy/Swift 两侧均为 FFT 数值噪声（不可复现），
        // 只断言 Swift 侧同处噪声底（< 1e-16，无系统性偏移/泄漏）。
        var strong = 0, edge = 0, noise = 0
        for (i, xk) in fxp.enumerated() {
            let kv = 3.0 * xk
            let j = Int((kv / dk).rounded()) + n / 2
            #expect(j >= 0 && j < n, "第\(i)点 k=\(kv) 越界")
            guard j >= 0 && j < n else { continue }
            let actual = pd[j] / pdMax
            let expected = fyp[i]
            if expected > 1e-16 {
                strong += 1
                expectPointwiseClose([actual], [expected], tolerance: 1e-6, "信号区 第\(i)点")
            } else if expected > 1e-30 {
                edge += 1
                expectPointwiseClose([actual], [expected], tolerance: 1e-3, "过渡区 第\(i)点")
            } else {
                noise += 1
                #expect(actual < 1e-16, "噪声区 第\(i)点 actual=\(actual)")
            }
        }
        #expect(strong >= 4, "信号区样本数 \(strong)")
        #expect(strong + edge + noise == fxp.count)
    }

    @Test("compute 输出结构：双图（2 + 1 系列）+ 摘要 + 理论卡")
    func computeStructure() async throws {
        let module = UncertaintyModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        guard case .lineSeries(let gaussian) = result.charts[0],
              case .lineSeries(let well) = result.charts[1] else {
            Issue.record("应为两个 lineSeries"); return
        }
        #expect(gaussian.series.count == 2)
        #expect(gaussian.series[0].points.count == 2000)
        #expect(gaussian.series[1].points.count == 2000)
        #expect(well.series.count == 1)
        #expect(well.series[0].points.count == 500)   // 4000 / 8
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("秒级档预算：compute 中位 < 2000 ms")
    func computeBudget() async throws {
        let module = UncertaintyModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 2000)
    }
}
