import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W5-2 全同粒子（费米气体 + HOM）：脚本打印值断言 + fixtures 曲线对拍。
@Suite("W5 全同粒子：费米气体与 HOM")
struct IdenticalParticlesModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    /// 脚本默认：Cu，n = 8.49e28 m⁻³。
    private func cu() throws -> (kF: Double, eF: Double, eFJ: Double) {
        let cs = try constants()
        let hbar = try cs.value("hbar"), me = try cs.value("m_e"), eV = try cs.value("eV")
        let kF = IdenticalParticlesMath.kF(8.49e28)
        return (kF, IdenticalParticlesMath.eF(8.49e28, hbar: hbar, me: me) / eV,
                IdenticalParticlesMath.eF(8.49e28, hbar: hbar, me: me))
    }

    @Test("铜的费米气体参数：kF/EF/TF/vF（脚本打印值，rel < 1e-3）")
    func fermiGasCu() throws {
        let cs = try constants()
        let kB = try cs.value("kB"), hbar = try cs.value("hbar"), me = try cs.value("m_e")
        let (kF, eFeV, eFJ) = try cu()
        #expect(abs(kF / 1.3597e10 - 1) < 1e-3, "k_F = \(kF)")
        #expect(abs(eFeV / 7.0438 - 1) < 1e-3, "E_F = \(eFeV) eV")
        #expect(abs(eFJ / kB / 8.1740e4 - 1) < 1e-3, "T_F")
        #expect(abs(hbar * kF / me / 1.5741e6 - 1) < 1e-3, "v_F")
    }

    @Test("态密度标度：D(E) ∝ E^(1/2)——4E 处 2 倍")
    func dosScaling() throws {
        let cs = try constants()
        let hbar = try cs.value("hbar"), me = try cs.value("m_e")
        let d1 = IdenticalParticlesMath.dos(1, hbar: hbar, me: me)
        let d4 = IdenticalParticlesMath.dos(4, hbar: hbar, me: me)
        #expect(abs(d4 / d1 - 2) < 1e-12)
    }

    @Test("HOM：τ=0 塌缩到 0，|τ|→∞ 趋 0.5；可见度 1 vs 0")
    func homDip() {
        let sigma = 1.0e12
        #expect(IdenticalParticlesMath.homRate(0, sigmaW: sigma) < 1e-20)
        #expect(abs(IdenticalParticlesMath.homRate(1e-6, sigmaW: sigma) - 0.5) < 1e-9)
        let tau = Num.linspace(-4e-12, 4e-12, count: 800)
        let rc = tau.map { IdenticalParticlesMath.homRate($0, sigmaW: sigma) }
        #expect(abs(IdenticalParticlesMath.visibility(rc) - 0.99995) < 1e-4,
                "V = \(IdenticalParticlesMath.visibility(rc))")
        #expect(IdenticalParticlesMath.visibility([Double](repeating: 0.5, count: 800)) == 0)
    }

    @Test("fixture 曲线对拍：D(E)（600 点网格口径）+ HOM 双曲线（800 点）")
    func curvesMatchFixtures() throws {
        let fx = try ModuleFixture.load("全同粒子_费米气体与HOM__v2022")
        let cs = try constants()
        let hbar = try cs.value("hbar"), me = try cs.value("m_e"), eV = try cs.value("eV")

        // fig0：D(E)，x = E/eV（脚本网格 linspace(1e-3, 3EF, 600) 的降采样）
        let (x, y) = try fx.line(0, 0, index: 0)
        #expect(x.count == 512)
        // 脚本网格口径：E = linspace(1e-3, 3·EF, 600)（J；上限为 3 倍费米能）
        let eF = IdenticalParticlesMath.eF(8.49e28, hbar: hbar, me: me)
        let eMax = 3 * eF
        for xi in x {
            let eJ = xi * eV
            #expect(eJ <= 1e-3 + 1e-12 && eJ >= eMax - 1e-18,
                    "x 应落在脚本网格 [1e-3, 3EF] 内（J 口径）")
        }
        expectPointwiseClose(x.map { IdenticalParticlesMath.dos($0 * eV, hbar: hbar, me: me) * eV },
                             y, tolerance: 1e-8, "D(E)∝E^1/2")
        // E_F 竖线（fixture line 1）
        let (xf, _) = try fx.line(0, 0, index: 1)
        #expect(abs(xf[0] - eF / eV) < 1e-8)

        // fig1：HOM 双曲线（x = τ ps）
        let (x1, y1) = try fx.line(1, 0, index: 0)
        let (_, y2) = try fx.line(1, 0, index: 1)
        #expect(x1.count == 512)
        expectPointwiseClose(x1.map { IdenticalParticlesMath.homRate($0 * 1e-12, sigmaW: 1e12) },
                             y1, tolerance: 1e-8, "全同光子 R_c(τ)")
        expectPointwiseClose([Double](repeating: 0.5, count: x1.count),
                             y2, tolerance: 1e-12, "非全同光子 R_c = 0.5")
    }

    @Test("compute 输出结构：2 图（1+2 系列）+ E_F 参考线 + 摘要 + 理论卡")
    func computeStructure() async throws {
        let module = IdenticalParticlesModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        guard case .lineSeries(let c0) = result.charts[0],
              case .lineSeries(let c1) = result.charts[1] else {
            Issue.record("应为 lineSeries"); return
        }
        #expect(c0.series.count == 1)
        #expect(c0.series[0].points.count == 300)   // 600 / 2
        #expect(c0.referenceLines.count == 1)
        #expect(c1.series.count == 2)
        #expect(c1.series[0].points.count == 400)   // 800 / 2
        #expect(result.summary.count == 5)
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("实时档预算：compute 中位 < 16 ms")
    func computeBudget() async throws {
        let module = IdenticalParticlesModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }
}
