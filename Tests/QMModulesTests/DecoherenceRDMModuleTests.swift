import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W5-1 退相干约化密度矩阵：fixtures 曲线对拍（300 点 ×3 曲线）+ 端点断言。
@Suite("W5 退相干约化密度矩阵")
struct DecoherenceRDMModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    @Test("端点：t=0 相干/纯度/熵；t→∞ 纯度 0.5、熵 1 ebit（脚本打印值）")
    func endpoints() {
        #expect(abs(DecoherenceRDMMath.offDiagonal(0, gamma: 1) - 1) < 1e-12)
        #expect(abs(DecoherenceRDMMath.purity(0, gamma: 1) - 1) < 1e-12)
        #expect(DecoherenceRDMMath.entropyEbits(0, gamma: 1) < 1e-9)
        #expect(abs(DecoherenceRDMMath.purity(100, gamma: 1) - 0.5) < 1e-12)
        #expect(abs(DecoherenceRDMMath.entropyEbits(100, gamma: 1) - 1) < 1e-9)
    }

    @Test("纯度与熵的单调性：Γt 增大 → 纯度降、熵升")
    func monotonicity() {
        var lastP = 2.0, lastS = -1.0
        for i in 0...40 {
            let t = Double(i) * 0.1
            let p = DecoherenceRDMMath.purity(t, gamma: 1)
            let s = DecoherenceRDMMath.entropyEbits(t, gamma: 1)
            #expect(p <= lastP + 1e-12)
            #expect(s >= lastS - 1e-12)
            lastP = p; lastS = s
        }
    }

    @Test("相干半衰期：e^(−Γt) = 0.5 ⟹ t = ln2/Γ")
    func halfLife() {
        let tHalf = log(2.0)
        #expect(abs(DecoherenceRDMMath.offDiagonal(tHalf, gamma: 1) - 0.5) < 1e-12)
    }

    @Test("fixture 曲线对拍：非对角元 / 纯度 / 熵（Γ=1，t ∈ [0,4] 300 点）")
    func curvesMatchFixtures() throws {
        let fx = try ModuleFixture.load("绘景变换与密度矩阵_退相干约化密度矩阵__v2022")
        let (xt, yOff) = try fx.line(0, 0, index: 0)
        let (_, yPur) = try fx.line(0, 0, index: 1)
        let (_, yEnt) = try fx.line(0, 1, index: 0)
        #expect(xt.count == 300 && yOff.count == 300)

        // 网格口径：fixture 存 9 位有效数字（舍入 ~5e-9），验证等差 + 端点
        #expect(abs(xt.first! - 0) < 1e-12 && abs(xt.last! - 4) < 1e-9)
        let step = 4.0 / 299
        for i in 1..<xt.count {
            #expect(abs(xt[i] - xt[i - 1] - step) < 5e-8, "网格 t 应为 linspace(0,4,300)")
        }
        expectPointwiseClose(xt.map { DecoherenceRDMMath.offDiagonal($0, gamma: 1) },
                             yOff, tolerance: 1e-8, "ρ₁₂(t)")
        expectPointwiseClose(xt.map { DecoherenceRDMMath.purity($0, gamma: 1) },
                             yPur, tolerance: 1e-8, "Tr ρ²")
        expectPointwiseClose(xt.map { DecoherenceRDMMath.entropyEbits($0, gamma: 1) },
                             yEnt, tolerance: 1e-8, "S(ρ_A)/ln2")
    }

    @Test("compute 输出结构：2 图（2+1 系列）+ 摘要 + 理论卡")
    func computeStructure() async throws {
        let module = DecoherenceRDMModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        guard case .lineSeries(let c0) = result.charts[0],
              case .lineSeries(let c1) = result.charts[1] else {
            Issue.record("应为 lineSeries"); return
        }
        #expect(c0.series.count == 2)
        #expect(c0.series[0].points.count == 150)   // 300 / 2
        #expect(c1.series.count == 1)
        #expect(c1.referenceLines.count == 1)
        #expect(result.summary.count == 3)
        #expect(result.theory?.formulas.count == 3)
    }

    @Test("实时档预算：compute 中位 < 16 ms")
    func computeBudget() async throws {
        let module = DecoherenceRDMModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }
}
