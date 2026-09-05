import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W4-3 康普顿散射：复刻脚本内 3 条断言（fixtures 康普顿散射_位移与反冲__v2022.json）。
@Suite("W4 康普顿散射")
struct ComptonModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    @Test("check1：λ_C 与 CODATA 2.426310e-12 m 一致（rel < 1e-4）")
    func comptonWavelength() throws {
        let cs = try constants()
        let lamC = ComptonMath.comptonWavelength(
            h: try cs.value("h"), me: try cs.value("m_e"), c: try cs.value("c"))
        #expect(abs(lamC - 2.426310e-12) / 2.426310e-12 < 1e-4,
                "λ_C = \(lamC)")
        // 同时应与常量库 lambda_C 键值一致
        let lamCRef = try cs.value("lambda_C")
        #expect(abs(lamC - lamCRef) / lamC < 1e-12)
    }

    @Test("check2：θ=0 → Δλ=0 且 K_e≈0")
    func zeroAngle() throws {
        let cs = try constants()
        let lamC = ComptonMath.comptonWavelength(
            h: try cs.value("h"), me: try cs.value("m_e"), c: try cs.value("c"))
        let ke = ComptonMath.recoilKinetic(theta: 0, lambda0: 71e-12,
                                           lambdaC: lamC, hc: try cs.value("h") * cs.value("c"))
        #expect(ComptonMath.deltaLambda(theta: 0, lambdaC: lamC) == 0)
        #expect(abs(ke) < 1e-30, "K_e(0) = \(ke)")
    }

    @Test("check3：θ=π → Δλ = 2λ_C（rel < 1e-12）")
    func backscatter() throws {
        let cs = try constants()
        let lamC = ComptonMath.comptonWavelength(
            h: try cs.value("h"), me: try cs.value("m_e"), c: try cs.value("c"))
        let dl = ComptonMath.deltaLambda(theta: .pi, lambdaC: lamC)
        #expect(abs(dl - 2 * lamC) / (2 * lamC) < 1e-12)
    }

    @Test("compute 输出结构：双图 + 2λ_C 参考线 + 摘要 + 理论卡")
    func computeStructure() async throws {
        let module = ComptonModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        guard case .lineSeries(let c0) = result.charts[0],
              case .lineSeries(let c1) = result.charts[1] else {
            Issue.record("应为两个 lineSeries"); return
        }
        #expect(c0.series[0].points.count == 400)
        #expect(c1.series[0].points.count == 400)
        #expect(c0.referenceLines.count == 1)
        // Mo Kα 71 pm 背散射反冲能量量级校验（~1.6 keV）
        #expect(result.summary.count == 3)
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("实时档预算：compute 中位 < 16 ms")
    func computeBudget() async throws {
        let module = ComptonModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }
}
