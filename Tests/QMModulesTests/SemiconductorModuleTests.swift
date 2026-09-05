import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W5-5 半导体：n_i(T) 与 Shockley I-V（脚本打印值断言；主曲线无 label
/// 未入 fixtures，公式复算对拍——W4 双缝同款口径）。
@Suite("W5 半导体：载流子与 pn 结")
struct SemiconductorModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    /// 脚本默认：E_g = 1.12 eV（Si）。
    private func egJ() throws -> Double {
        try constants().value("eV") * 1.12
    }

    @Test("n_i(T)：脚本网格点三档（T = 300.50/399.50/600.00 K，rel < 1e-3）")
    func niValues() throws {
        let kB = try constants().value("kB")
        let eg = try egJ()
        // 脚本网格 T = linspace(100, 600, 400)，打印取最近点：
        // idx 160 → 300.5013 / idx 239 → 399.4987 / idx 399 → 600.0
        let cases: [(Double, Double)] = [(300.501253, 6.9389e9), (399.498747, 2.2602e12),
                                         (600.0, 9.5466e14)]
        for (T, exp) in cases {
            let ni = SemiconductorMath.ni(T, eg: eg, kB: kB)
            #expect(abs(ni / exp - 1) < 1e-3, "n_i(\(T)K) = \(ni)")
        }
    }

    @Test("n_i 指数主导：E_g 每 +0.1 eV，n_i(300 K) 掉接近一个量级")
    func niEgDominance() throws {
        let kB = try constants().value("kB")
        let eV = try constants().value("eV")
        let ni1 = SemiconductorMath.ni(300, eg: 1.12 * eV, kB: kB)
        let ni2 = SemiconductorMath.ni(300, eg: 1.22 * eV, kB: kB)
        // exp(0.1 eV / 2k_BT) ≈ exp(1.93) ≈ 6.9 倍
        #expect(ni1 / ni2 > 6, "ratio = \(ni1 / ni2)")
        #expect(ni1 / ni2 < 8, "ratio = \(ni1 / ni2)")
    }

    @Test("热电压 k_BT/e = 25.85 mV @300 K；二极管电流四档 + 反向饱和")
    func diodeIV() throws {
        let cs = try constants()
        let kB = try cs.value("kB"), e = try cs.value("e")
        let vT = SemiconductorMath.thermalVoltage(T: 300, kB: kB, e: e)
        #expect(abs(vT * 1000 - 25.85) < 0.01, "\(vT * 1000) mV")

        let isat = 1e-12
        let expected: [(Double, Double)] = [(0.5, 2.510e-4), (0.6, 1.201e-2),
                                             (0.7, 5.748e-1), (0.8, 2.750e1)]
        for (v, iExp) in expected {
            let i = SemiconductorMath.diodeCurrent(v, isat: isat, vT: vT)
            #expect(abs(i / iExp - 1) < 1e-3, "I(\(v)V) = \(i)")
        }
        let iRev = SemiconductorMath.diodeCurrent(-1.0, isat: isat, vT: vT)
        #expect(abs(iRev + 1e-12) < 1e-15, "反向饱和 −I_s")
    }

    @Test("I-V 指数律：正向每 60 mV（≈ vT·ln2×2... 实测 ln10·vT）电流 ×10")
    func decadeSlope() throws {
        let cs = try constants()
        let kB = try cs.value("kB"), e = try cs.value("e")
        let vT = SemiconductorMath.thermalVoltage(T: 300, kB: kB, e: e)
        let isat = 1e-12
        let i1 = SemiconductorMath.diodeCurrent(0.6, isat: isat, vT: vT)
        let i2 = SemiconductorMath.diodeCurrent(0.6 + vT * log(10), isat: isat, vT: vT)
        #expect(abs(i2 / i1 / 10 - 1) < 1e-9, "ln10·vT ≈ 59.6 mV/十倍程")
    }

    @Test("compute 输出结构：2 图（400 点 + 450 点）+ 参考线 + 摘要 + 理论卡")
    func computeStructure() async throws {
        let module = SemiconductorModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        guard case .lineSeries(let c0) = result.charts[0],
              case .lineSeries(let c1) = result.charts[1] else {
            Issue.record("应为 lineSeries"); return
        }
        #expect(c0.series[0].points.count == 400)
        #expect(c0.spec.yAxis.scale == .log)
        #expect(c1.series[0].points.count == 450)   // 900 / 2
        #expect(c1.spec.yAxis.scale == .log)
        #expect(c1.referenceLines.count == 2)
        #expect(result.summary.count == 3)
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("实时档预算：compute 中位 < 16 ms")
    func computeBudget() async throws {
        let module = SemiconductorModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }
}
