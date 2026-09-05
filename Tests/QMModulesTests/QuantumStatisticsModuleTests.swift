import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W5-6 量子统计：BEC T_c/凝聚分数（脚本打印值）+ 费米简并压
/// + fixtures P₀(n) 曲线对拍。
@Suite("W5 量子统计：BEC 与费米气体")
struct QuantumStatisticsModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    @Test("Rb-87 BEC：T_c = 8.58e-8 K（n = 1e19，脚本打印值，rel < 1e-3）")
    func becTc() throws {
        let cs = try constants()
        let hbar = try cs.value("hbar"), kB = try cs.value("kB"), amu = try cs.value("u")
        let tc = QuantumStatisticsMath.becTc(n: 1e19, m: 86.909 * amu, hbar: hbar, kB: kB)
        #expect(abs(tc / 8.58e-8 - 1) < 1e-3, "T_c = \(tc) K")
        // 标度律：m 翻倍 → T_c × 1/2；n 翻倍 → T_c × 2^(2/3)
        let tc2 = QuantumStatisticsMath.becTc(n: 1e19, m: 2 * 86.909 * amu, hbar: hbar, kB: kB)
        #expect(abs(tc2 / tc - 0.5) < 1e-12)
        let tc3 = QuantumStatisticsMath.becTc(n: 2e19, m: 86.909 * amu, hbar: hbar, kB: kB)
        #expect(abs(tc3 / tc - pow(2, 2.0 / 3)) < 1e-12)
    }

    @Test("凝聚分数：1−(T/T_c)^1.5；T_c 以上为 0；t=0.5 时 1−2^(−3/2)")
    func condensateFraction() {
        #expect(abs(QuantumStatisticsMath.condensateFraction(0) - 1) < 1e-12)
        #expect(abs(QuantumStatisticsMath.condensateFraction(0.5) - (1 - pow(0.5, 1.5))) < 1e-12)
        #expect(QuantumStatisticsMath.condensateFraction(1) == 0)
        #expect(QuantumStatisticsMath.condensateFraction(1.3) == 0)
    }

    @Test("铜：E_F = 7.044 eV，P₀ = 3.833e10 Pa（脚本打印值）")
    func copper() throws {
        let cs = try constants()
        let hbar = try cs.value("hbar"), me = try cs.value("m_e"), eV = try cs.value("eV")
        let ef = QuantumStatisticsMath.fermiEnergy(8.49e28, hbar: hbar, me: me)
        let p0 = QuantumStatisticsMath.degeneracyPressure(8.49e28, eF: ef)
        #expect(abs(ef / eV / 7.044 - 1) < 1e-3, "E_F = \(ef / eV) eV")
        #expect(abs(p0 / 3.833e10 - 1) < 1e-3, "P₀ = \(p0) Pa")
    }

    @Test("白矮星：E_F ≈ 164 keV，P₀ = 3.161e21 Pa（脚本打印值）")
    func whiteDwarf() throws {
        let cs = try constants()
        let hbar = try cs.value("hbar"), me = try cs.value("m_e"), eV = try cs.value("eV"), amu = try cs.value("u")
        let nWd = 1e9 / (2 * amu)
        let ef = QuantumStatisticsMath.fermiEnergy(nWd, hbar: hbar, me: me)
        let p0 = QuantumStatisticsMath.degeneracyPressure(nWd, eF: ef)
        #expect(abs(ef / eV / 163813.95 - 1) < 1e-3, "E_F = \(ef / eV) eV")
        #expect(abs(p0 / 3.161e21 - 1) < 1e-3, "P₀ = \(p0) Pa")
    }

    @Test("P₀ ∝ n^(5/3)：loglog 斜率 = 5/3")
    func pressureSlope() throws {
        let cs = try constants()
        let hbar = try cs.value("hbar"), me = try cs.value("m_e")
        let n1 = 1e28, n2 = 1e35
        let p1 = QuantumStatisticsMath.degeneracyPressure(
            n1, eF: QuantumStatisticsMath.fermiEnergy(n1, hbar: hbar, me: me))
        let p2 = QuantumStatisticsMath.degeneracyPressure(
            n2, eF: QuantumStatisticsMath.fermiEnergy(n2, hbar: hbar, me: me))
        let slope = log(p2 / p1) / log(n2 / n1)
        #expect(abs(slope - 5.0 / 3) < 1e-9, "slope = \(slope)")
    }

    @Test("fixture 曲线对拍：P₀(n)（logspace(27,36) 400 点）")
    func curvesMatchFixtures() throws {
        let fx = try ModuleFixture.load("量子统计_BEC与费米气体__v2022")
        let cs = try constants()
        let hbar = try cs.value("hbar"), me = try cs.value("m_e")
        let (x, y) = try fx.line(1, 0, index: 0)
        #expect(x.count == 400)
        // 脚本网格：n = logspace(27, 36, 400)（fixture 存 9 位有效数字，相对舍入 ~5e-9）
        for i in x.indices {
            #expect(abs(x[i] - pow(10, 27 + 9.0 * Double(i) / 399)) / x[i] < 1e-8,
                    "n 网格应为 logspace(27,36,400)")
        }
        let p = x.map {
            QuantumStatisticsMath.degeneracyPressure(
                $0, eF: QuantumStatisticsMath.fermiEnergy($0, hbar: hbar, me: me))
        }
        expectPointwiseClose(p, y, tolerance: 1e-8, "P₀ ∝ n^(5/3)")
    }

    @Test("compute 输出结构：2 图（400 点 ×2）+ 摘要 + 理论卡")
    func computeStructure() async throws {
        let module = QuantumStatisticsModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        guard case .lineSeries(let c0) = result.charts[0],
              case .lineSeries(let c1) = result.charts[1] else {
            Issue.record("应为 lineSeries"); return
        }
        #expect(c0.series[0].points.count == 400)
        #expect(c1.series[0].points.count == 400)
        #expect(c0.referenceLines.count == 1)
        #expect(c1.referenceLines.count == 2)   // 铜 + 白矮星
        #expect(result.summary.count == 3)
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("实时档预算：compute 中位 < 16 ms")
    func computeBudget() async throws {
        let module = QuantumStatisticsModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }
}
