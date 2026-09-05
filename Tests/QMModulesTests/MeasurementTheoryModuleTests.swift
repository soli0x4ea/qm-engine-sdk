import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W4-7 指针耦合与芝诺：复刻脚本 4 条断言 + fixtures 曲线对拍（np.interp 口径复刻）。
@Suite("W4 指针耦合与芝诺")
struct MeasurementTheoryModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }
    private static let sigma = 1.0
    private static let alpha2 = 0.6
    private static let lambdas = [0.0, 1.0, 2.0, 4.0]

    private func grid() -> (x: [Double], phi0: [Double]) {
        MeasurementTheoryMath.pointerGrid(sigma: Self.sigma)
    }

    @Test("check1：指针初态归一化 < 1e-6（trapz 复刻）")
    func pointerNormalized() {
        let (x, phi0) = grid()
        let dx = x[1] - x[0]
        var sq = [Double](repeating: 0, count: x.count)
        for i in x.indices { sq[i] = phi0[i] * phi0[i] }
        let nrm = MeasurementTheoryMath.trapz(sq, dx: dx)
        #expect(abs(nrm - 1.0) < 1e-6, "nrm = \(nrm)")
    }

    @Test("check2：重叠积分数值 vs 解析 e^(−λ²) 最大偏差 < 2e-2")
    func overlapDecay() {
        let (x, phi0) = grid()
        for lam in Self.lambdas {
            let (_, overlap) = MeasurementTheoryMath.jointProbability(
                lambda: lam, sigma: Self.sigma, alpha2: Self.alpha2, grid: x, phi0v: phi0)
            let ana = exp(-lam * lam)
            #expect(abs(overlap - ana) < 2e-2,
                    "λ=\(lam): num=\(overlap) ana=\(ana)")
        }
    }

    @Test("check3/4：芝诺 P(N=1) < 1e-9；P(N=512) > 0.99")
    func zenoBounds() {
        let p1 = MeasurementTheoryMath.zenoSurvival(N: 1, omega: 1.0, tTotal: .pi)
        let p512 = MeasurementTheoryMath.zenoSurvival(N: 512, omega: 1.0, tTotal: .pi)
        #expect(abs(p1) < 1e-9, "P(1) = \(p1)")
        #expect(p512 > 0.99, "P(512) = \(p512)")
    }

    @Test("fixture 曲线对拍：4 条 λ 曲线 512 点（interp 复刻，rel < 1e-6）")
    func curvesMatchFixtures() throws {
        let fx = try ModuleFixture.load("测量理论_指针耦合与芝诺__v2022")
        let (x, phi0) = grid()
        for (i, lam) in Self.lambdas.enumerated() {
            let (fx_, fy) = try fx.line(0, 0, index: i)
            // 在 fixture 的 x 网格上按 np.interp 口径重算（脚本 y = prob/max(prob)，max 取 4000 网格）
            var plus = [Double](repeating: 0, count: fx_.count)
            var minus = [Double](repeating: 0, count: fx_.count)
            for j in fx_.indices {
                plus[j] = MeasurementTheoryMath.interpShiftedPhi(fx_[j], grid: x, phi0v: phi0, lambda: lam)
                minus[j] = MeasurementTheoryMath.interpShiftedPhi(fx_[j], grid: x, phi0v: phi0, lambda: -lam)
            }
            var fullPlus = [Double](repeating: 0, count: x.count)
            var fullMinus = [Double](repeating: 0, count: x.count)
            for j in x.indices {
                fullPlus[j] = MeasurementTheoryMath.interpShiftedPhi(x[j], grid: x, phi0v: phi0, lambda: lam)
                fullMinus[j] = MeasurementTheoryMath.interpShiftedPhi(x[j], grid: x, phi0v: phi0, lambda: -lam)
            }
            var peak = 0.0
            var probs = [Double](repeating: 0, count: x.count)
            for j in x.indices {
                probs[j] = Self.alpha2 * fullPlus[j] * fullPlus[j]
                    + (1 - Self.alpha2) * fullMinus[j] * fullMinus[j]
                if probs[j] > peak { peak = probs[j] }
            }
            let actual = (0..<fx_.count).map { j in
                (Self.alpha2 * plus[j] * plus[j] + (1 - Self.alpha2) * minus[j] * minus[j]) / peak
            }
            expectPointwiseClose(actual, fy, tolerance: 1e-6, "指针耦合 λ=\(lam)")
        }
    }

    @Test("compute 输出结构：双图（4 系列 + semilogx 10 点）+ 摘要 + 理论卡")
    func computeStructure() async throws {
        let module = MeasurementTheoryModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        guard case .lineSeries(let pointer) = result.charts[0],
              case .lineSeries(let zeno) = result.charts[1] else {
            Issue.record("应为两个 lineSeries"); return
        }
        #expect(pointer.series.count == 4)
        #expect(pointer.series[0].points.count == 2000)   // 4000 / 2
        #expect(zeno.series[0].points.count == 10)
        #expect(zeno.spec.xAxis.scale == .log)             // semilogx 首用
        #expect(zeno.referenceLines.count == 1)
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("芝诺单调性：N 增大 → P(N) 单调不减（Ω=1）")
    func zenoMonotone() {
        let Ns = MeasurementTheoryMath.zenoCounts
        let Ps = Ns.map { MeasurementTheoryMath.zenoSurvival(N: $0, omega: 1.0, tTotal: .pi) }
        for i in 1..<Ps.count {
            #expect(Ps[i] >= Ps[i - 1] - 1e-15, "N=\(Ns[i]) P=\(Ps[i])")
        }
    }

    /// naive 参考实现：收尾包 1 快路径**前**的 jointProbability 旧路径逐句复刻。
    /// 用于锁死 λ=0 复用与不变量外提的 bit 级无损声明。
    private static func jointProbabilityLegacy(
        lambda: Double, sigma: Double, alpha2: Double, grid: [Double], phi0v: [Double]
    ) -> (y: [Double], overlap: Double) {
        let n = grid.count
        let dx = grid[1] - grid[0]
        let beta2 = 1 - alpha2
        var prob = [Double](repeating: 0, count: n)
        var sum = 0.0
        var peak = 0.0
        for i in 0..<n {
            let xq = grid[i]
            let plus = MeasurementTheoryMath.interpShiftedPhi(xq, grid: grid, phi0v: phi0v, lambda: lambda)
            let minus = MeasurementTheoryMath.interpShiftedPhi(xq, grid: grid, phi0v: phi0v, lambda: -lambda)
            sum += plus * minus
            let p = alpha2 * plus * plus + beta2 * minus * minus
            prob[i] = p
            if p > peak { peak = p }
        }
        let first = MeasurementTheoryMath.interpShiftedPhi(grid[0], grid: grid, phi0v: phi0v, lambda: lambda)
            * MeasurementTheoryMath.interpShiftedPhi(grid[0], grid: grid, phi0v: phi0v, lambda: -lambda)
        let last = MeasurementTheoryMath.interpShiftedPhi(grid[n - 1], grid: grid, phi0v: phi0v, lambda: lambda)
            * MeasurementTheoryMath.interpShiftedPhi(grid[n - 1], grid: grid, phi0v: phi0v, lambda: -lambda)
        let overlap = dx * (sum - 0.5 * (first + last))
        for i in 0..<n { prob[i] /= peak }
        return (prob, overlap)
    }

    @Test("jointProbability 快路径 bit 级无损：λ=0 与 λ≠0 均与旧路径逐元素 ==（收尾包 1）")
    func jointProbabilityBitwiseStable() {
        let (x, phi0) = grid()
        for lam in [0.0, 1.0, 2.0, 4.0] {
            let legacy = Self.jointProbabilityLegacy(
                lambda: lam, sigma: Self.sigma, alpha2: Self.alpha2, grid: x, phi0v: phi0)
            let current = MeasurementTheoryMath.jointProbability(
                lambda: lam, sigma: Self.sigma, alpha2: Self.alpha2, grid: x, phi0v: phi0)
            #expect(legacy.y.count == current.y.count)
            for i in x.indices {
                #expect(legacy.y[i] == current.y[i],
                        "λ=\(lam) y[\(i)]: \(legacy.y[i]) vs \(current.y[i])")
            }
            #expect(legacy.overlap == current.overlap,
                    "λ=\(lam) overlap: \(legacy.overlap) vs \(current.overlap)")
        }
    }

    @Test("实时档预算：compute 中位 < 16 ms")
    func computeBudget() async throws {
        let module = MeasurementTheoryModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }
}
