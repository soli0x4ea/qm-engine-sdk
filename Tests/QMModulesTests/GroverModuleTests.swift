import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W6 笔记 34《量子计算与量子算法》：Grover 搜索成功概率。
/// fixtures：量子计算_Grover搜索__v2022.json（P(n) 61 点 N=1024 + 最优迭代 n*=25 竖线）。
@Suite("W6 Grover 搜索")
struct GroverModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    private static let N = 1024.0

    // MARK: - 物理律

    @Test("物理律：P(0) = 1/N（均匀叠加起点）；sin² 律精确成立")
    func initialProbabilityLaw() {
        #expect(abs(GroverMath.probability(n: 0, N: Self.N) - 1.0 / Self.N) < 1e-15,
                "n=0 时 P = sin²θ = 1/N = 1/1024")
        // sin² 律：P(n) = sin²((2n+1)θ)，θ = arcsin(1/√N)
        let theta = GroverMath.theta(N: Self.N)
        for n in [0, 1, 7, 25, 60] {
            let a = sin((2.0 * Double(n) + 1.0) * theta)
            #expect(abs(GroverMath.probability(n: n, N: Self.N) - a * a) < 1e-15)
        }
        #expect(abs(GroverMath.theta(N: Self.N) - asin(1.0 / 32.0)) < 1e-15)
    }

    @Test("物理律：峰值恰在 n* = 25（sin² 律 + √N 预测）；P(n*) ≈ 1；过冲回落")
    func peakLaw() {
        let probs = GroverMath.probabilities(N: Self.N, nMax: 60)
        #expect(probs.count == 61)
        // 峰值位置 = 网格 argmax = n* = 25
        let argmax = probs.indices.max { probs[$0] < probs[$1] } ?? -1
        #expect(argmax == 25, "argmax = \(argmax)")
        #expect(GroverMath.optimalIterations(N: Self.N) == 25)
        #expect(GroverMath.peakProbability(N: Self.N) == probs[25])
        #expect(probs[25] > 0.999, "峰值概率 P(25) = \(probs[25])")
        // √N 预测：n* ≈ (π/4)√N = 25.13，与整数 n* 距离 < 1
        #expect(abs(Double(argmax) - Double.pi / 4.0 * 32.0) < 1.0)
        // sin² 律的攀升 / 过冲：n < n* 单调升；n* < n < π/(2θ)−1/2 回落至低谷
        for n in 0..<25 { #expect(probs[n + 1] > probs[n], "攀升段 n=\(n)") }
        for n in 25..<50 { #expect(probs[n + 1] < probs[n], "过冲回落 n=\(n)") }
        // 低谷在 (2n+1)θ = π 即 n ≈ 49.7（网格 n=50），其后 sin² 周期性回升
        #expect(probs[51] > probs[50], "低谷 n=50 后周期回升")
    }

    @Test("物理律：查询复杂度——量子 (π/4)√N vs 经典 N/2，加速比 2√N/π")
    func speedupLaw() {
        #expect(abs(GroverMath.quantumQueries(N: Self.N) - Double.pi / 4.0 * 32.0) < 1e-12)
        #expect(abs(GroverMath.classicalQueries(N: Self.N) - 512) < 1e-12)
        #expect(abs(GroverMath.speedup(N: Self.N) - 2.0 * 32.0 / Double.pi) < 1e-12,
                "加速比 = 2√N/π ≈ 20.4")
        #expect(GroverMath.quantumQueries(N: Self.N) < GroverMath.classicalQueries(N: Self.N))
        // 加速比随 N 增长（平方根加速的本质）
        #expect(GroverMath.speedup(N: 1_048_576.0) > GroverMath.speedup(N: Self.N))
    }

    // MARK: - fixtures 对拍

    @Test("fixtures：P(n) 61 点逐点对拍（< 1e-8）+ 最优迭代竖线 n* = 25")
    func matchesFixtures() throws {
        let fx = try ModuleFixture.load("量子计算_Grover搜索__v2022")
        let (ns, pF) = try fx.line(0, 0, label: "sin^2")
        #expect(ns.count == 61)
        #expect(ns.map(Int.init) == Array(0...60))
        expectPointwiseClose(GroverMath.probabilities(N: Self.N, nMax: 60), pF,
                             tolerance: 1e-8, "P(n) = sin²((2n+1)θ)")

        let (xOpt, _) = try fx.line(0, 0, label: "optimal")
        #expect(xOpt[0] == 25 && xOpt[1] == 25, "最优迭代标注线 x = n* = 25")
        #expect(Double(xOpt[0]) == Double(GroverMath.optimalIterations(N: Self.N)))
    }

    // MARK: - compute 结构与预算

    @Test("compute 结构：1 图 2 系列（61 点）+ 峰值参考线 + 摘要 + 理论卡；实时档预算")
    func computeStructureAndBudget() async throws {
        let module = GroverModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 1)
        guard case .lineSeries(let c) = result.charts[0] else {
            Issue.record("应为 lineSeries"); return
        }
        #expect(c.series[0].points.count == 61)
        #expect(c.series[1].points.count == 1, "峰值标记点")
        #expect(c.referenceLines.count == 2)
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 4)

        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }
}
