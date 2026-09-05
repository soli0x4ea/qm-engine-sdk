import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W4-6 双缝干涉：复刻脚本计算（λ/dx 打印值）+ 干涉结构与周期性断言。
/// fixtures 曲线因脚本线未设 label 而未捕获（导出管线 _child 过滤），改为复刻脚本数值断言。
@Suite("W4 双缝干涉")
struct DoubleSlitModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    /// 脚本默认参数：E = 100 eV，d = 100 nm，L = 1 m。
    private func defaults() throws -> (lam: Double, fringe: Double) {
        let cs = try constants()
        let e = try cs.value("e"), h = try cs.value("h"), me = try cs.value("m_e")
        let lam = DoubleSlitMath.deBroglieWavelength(energyJoule: 100 * e, h: h, me: me)
        let fringe = DoubleSlitMath.fringeSpacing(lambda: lam, d: 1e-7, L: 1.0)
        return (lam, fringe)
    }

    @Test("λ(100 eV 电子) = 1.227 Å（脚本打印值，rel < 1e-3）")
    func deBroglie() throws {
        let (lam, _) = try defaults()
        #expect(abs(lam * 1e10 - 1.227) / 1.227 < 1e-3, "λ = \(lam * 1e10) Å")
    }

    @Test("条纹间距 Δx = λL/d = 1.227 mm（脚本打印值，rel < 1e-3）")
    func fringeSpacing() throws {
        let (_, fringe) = try defaults()
        #expect(abs(fringe * 1e3 - 1.2266) / 1.2266 < 1e-3, "Δx = \(fringe * 1e3) mm")
    }

    @Test("强度结构：中央极大 I(0) = 4I₀；周期性 I(Δx) = I(0)；半整数处归零")
    func intensityStructure() throws {
        let (lam, fringe) = try defaults()
        let d = 1e-7, L = 1.0
        let i0 = DoubleSlitMath.intensity(x: 0, lambda: lam, d: d, L: L)
        let iFringe = DoubleSlitMath.intensity(x: fringe, lambda: lam, d: d, L: L)
        let iHalf = DoubleSlitMath.intensity(x: fringe / 2, lambda: lam, d: d, L: L)
        #expect(abs(i0 - 4.0) < 1e-12, "I(0) = \(i0)")
        #expect(abs(iFringe - i0) < 1e-9, "周期性 I(Δx) = \(iFringe)")
        #expect(iHalf < 1e-12, "半整数条纹 I = \(iHalf)")
    }

    @Test("缝距反比：d 减半 → 条纹间距翻倍")
    func fringeInverseD() throws {
        let (lam, _) = try defaults()
        let f1 = DoubleSlitMath.fringeSpacing(lambda: lam, d: 1e-7, L: 1.0)
        let f2 = DoubleSlitMath.fringeSpacing(lambda: lam, d: 0.5e-7, L: 1.0)
        #expect(abs(f2 / f1 - 2) < 1e-12)
    }

    @Test("compute 输出结构：1 图 1 系列(500 点) + 2 参考线 + 摘要 + 理论卡")
    func computeStructure() async throws {
        let module = DoubleSlitModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 1)
        guard case .lineSeries(let c) = result.charts[0] else {
            Issue.record("应为 lineSeries"); return
        }
        #expect(c.series.count == 1)
        #expect(c.series[0].points.count == 500)   // 4000 / 8
        #expect(c.referenceLines.count == 2)
        // 参考线在 ±Δx 处
        let vals = c.referenceLines.map(\.value).sorted()
        #expect(abs(vals[0] + 1.2266) < 0.01 && abs(vals[1] - 1.2266) < 0.01)
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("实时档预算：compute 中位 < 16 ms")
    func computeBudget() async throws {
        let module = DoubleSlitModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }
}
