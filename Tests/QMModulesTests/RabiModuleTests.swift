import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W4-18 两能级 Rabi：复刻脚本（闭式 + RK4 绘景等价）+ fixtures 双线对拍。
@Suite("W4 两能级 Rabi")
struct RabiModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }
    private static let omega = 1.0
    private static let bigOmega = 2.0

    @Test("绘景等价：max|⟨σz⟩S − ⟨σz⟩I| < 1e-4（脚本打印口径）")
    func pictureEquivalence() {
        let (_, _, maxDiff) = RabiMath.sigmaZInteraction(
            omega: Self.omega, bigOmega: Self.bigOmega)
        #expect(maxDiff < 1e-4, "max diff = \(maxDiff)")
        #expect(maxDiff > 0, "RK4 与闭式应有非零但微小的步长误差")
    }

    @Test("RK4 曲线端点语义：t = 8π 处 ⟨σz⟩I ≈ 闭式（脚本 y[-1] = 0.951456）")
    func interactionEndpoint() {
        let (ts, expI, _) = RabiMath.sigmaZInteraction(
            omega: Self.omega, bigOmega: Self.bigOmega)
        let closed = RabiMath.sigmaZSchrodinger(
            t: ts[ts.count - 1], omega: Self.omega, bigOmega: Self.bigOmega)
        #expect(abs(expI[ts.count - 1] - closed) < 1e-5,
                "端点 \(expI[ts.count - 1]) vs 闭式 \(closed)")
        #expect(abs(expI[0] - 1) < 1e-15, "t=0 处 ⟨σz⟩ = 1")
    }

    @Test("Rabi 结构：Ω_R = √5；翻转深度 = (ω²−Ω²)/(ω²+Ω²) = −0.6；共振 ω=0 完全翻转")
    func rabiStructure() {
        let oR = RabiMath.rabiFrequency(omega: Self.omega, bigOmega: Self.bigOmega)
        #expect(abs(oR - 5.0.squareRoot()) < 1e-15)
        // 半周期 t = π/Ω_R 处 ⟨σz⟩ 取极值
        let tHalf = Double.pi / oR
        let z = RabiMath.sigmaZSchrodinger(t: tHalf, omega: Self.omega, bigOmega: Self.bigOmega)
        #expect(abs(z - (-0.6)) < 1e-12, "翻转深度 = \(z)")
        // 共振：ω = 0 ⇒ ⟨σz⟩ = cos²−sin² = cos(Ω t)，t = π/Ω 处到 −1
        let zRes = RabiMath.sigmaZSchrodinger(t: .pi / Self.bigOmega, omega: 0, bigOmega: Self.bigOmega)
        #expect(abs(zRes + 1) < 1e-12, "共振完全翻转 \(zRes)")
        // 周期性：⟨σz⟩(t + 2π/Ω_R) = ⟨σz⟩(t)
        let t1 = 0.7
        let z1 = RabiMath.sigmaZSchrodinger(t: t1, omega: Self.omega, bigOmega: Self.bigOmega)
        let z2 = RabiMath.sigmaZSchrodinger(t: t1 + 2 * .pi / oR, omega: Self.omega, bigOmega: Self.bigOmega)
        #expect(abs(z1 - z2) < 1e-12)
    }

    @Test("fixture 对拍：薛绘景闭式 400 点（rel < 1e-8）+ RK4 线（|差| < 1e-6）")
    func curvesMatchFixtures() throws {
        let fx = try ModuleFixture.load("绘景变换与密度矩阵_两能级Rabi__v2022")
        let (sxF, syF) = try fx.line(0, 0, label: "Schrödinger")
        let (_, iyF) = try fx.line(0, 0, label: "Interaction")
        #expect(sxF.count == 400 && iyF.count == 400)

        let (ts, expI, _) = RabiMath.sigmaZInteraction(
            omega: Self.omega, bigOmega: Self.bigOmega)
        for i in ts.indices {
            #expect(abs(sxF[i] - ts[i]) / (8 * .pi) < 1e-8, "t 网格 第\(i)点")
            let closed = RabiMath.sigmaZSchrodinger(
                t: ts[i], omega: Self.omega, bigOmega: Self.bigOmega)
            #expect(abs(syF[i] - closed) < 1e-8, "闭式 第\(i)点 \(syF[i]) vs \(closed)")
            #expect(abs(iyF[i] - expI[i]) < 1e-6, "RK4 第\(i)点 \(iyF[i]) vs \(expI[i])")
        }
    }

    @Test("compute 输出结构：1 图 2 系列（400 点）+ 摘要 + 理论卡")
    func computeStructure() async throws {
        let module = RabiModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 1)
        guard case .lineSeries(let c) = result.charts[0] else {
            Issue.record("应为 lineSeries"); return
        }
        #expect(c.series.count == 2)
        #expect(c.series[0].points.count == 400)
        #expect(c.series[1].points.count == 400)
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("秒级档预算：compute 中位 < 2000 ms")
    func computeBudget() async throws {
        let module = RabiModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 2000)
    }
}
