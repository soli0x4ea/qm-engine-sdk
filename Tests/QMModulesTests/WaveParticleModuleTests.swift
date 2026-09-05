import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W4-2 波粒二象性：fixtures 波粒二象性_波长与累积模拟__v2022.json 三条断言
/// + 模型 A 曲线逐点对拍（fixture 存 500 点 loglog 双曲线）。
@Suite("W4 波粒二象性")
struct WaveParticleModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    @Test("模型 A 曲线逐点对拍 fixtures（非相对论 + 相对论，500 点）")
    func curvesMatchFixtures() throws {
        let fx = try ModuleFixture.load("波粒二象性_波长与累积模拟__v2022")
        let cs = try constants()
        let h = try cs.value("h"), me = try cs.value("m_e")
        let e = try cs.value("e"), c = try cs.value("c")

        // fixture lines[0]=nonrelativistic, lines[1]=relativistic
        // （"nonrelativistic" 含子串 "relativistic"，标签匹配歧义，按下标取）
        let (nx, ny) = try fx.line(0, 0, index: 0)
        expectPointwiseClose(
            nx.map { WaveParticleMath.lambdaNonrel($0, h: h, me: me, e: e) * 1e9 },
            ny, tolerance: 1e-8, "非相对论 λ(V)")

        let (rx, ry) = try fx.line(0, 0, index: 1)
        expectPointwiseClose(
            rx.map { WaveParticleMath.lambdaRel($0, h: h, me: me, e: e, c: c) * 1e9 },
            ry, tolerance: 1e-8, "相对论 λ(V)")
    }

    @Test("脚本断言 1：54 V → 0.1669 nm（Davisson-Germer 主峰）")
    func davissonGermer() throws {
        let cs = try constants()
        let lam = WaveParticleMath.lambdaNonrel(
            54, h: try cs.value("h"), me: try cs.value("m_e"), e: try cs.value("e"))
        #expect(abs(lam * 1e9 - 0.1669) < 0.001, "λ(54V)=\(lam * 1e9) nm")
    }

    @Test("脚本断言 2：50 kV 相对论修正在 −3% ~ −1.5% 之间")
    func relativisticCorrection() throws {
        let cs = try constants()
        let h = try cs.value("h"), me = try cs.value("m_e")
        let e = try cs.value("e"), c = try cs.value("c")
        let corr = (WaveParticleMath.lambdaRel(50e3, h: h, me: me, e: e, c: c)
                    / WaveParticleMath.lambdaNonrel(50e3, h: h, me: me, e: e) - 1) * 100
        #expect(-3.0 < corr && corr < -1.5, "修正 \(corr)%")
    }

    @Test("脚本断言 3：N=60000 抽样与理论分布卡方一致（|z| < 3，种子 20260901）")
    func chiSquareConsistency() {
        let (x, p12) = WaveParticleMath.theoreticalP12()
        let hits = WaveParticleMath.sampleHits(x: x, p: p12, count: 60000, seed: 20260901)
        let r = WaveParticleMath.chiSquare(hits: hits, x: x, p: p12)
        #expect(abs(r.z) < 3, "z=\(r.z) (χ²=\(r.chi2), dof=\(r.dof))")
    }

    @Test("抽样可复现：同种子同序列；异种子异序列")
    func samplingReproducible() {
        let (x, p12) = WaveParticleMath.theoreticalP12()
        let a = WaveParticleMath.sampleHits(x: x, p: p12, count: 100, seed: 20260901)
        let b = WaveParticleMath.sampleHits(x: x, p: p12, count: 100, seed: 20260901)
        let c = WaveParticleMath.sampleHits(x: x, p: p12, count: 100, seed: 42)
        #expect(a == b)
        #expect(a != c)
    }

    @Test("实时档预算：N=60000 时 compute 中位 < 16 ms")
    func computeBudget() async throws {
        let module = WaveParticleModule()
        var values = ParamValues.defaults(for: module.params)
        values.discretes["N"] = "60000"
        try await expectComputeUnderBudget(
            module: module, values: values,
            constants: try constants(), budgetMillis: 16)
    }
}
