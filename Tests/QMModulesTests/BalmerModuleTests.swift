import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W4-4 巴尔末系：复刻脚本断言（R_H 精度 + 4 条谱线对观测值 rel < 1e-3）。
@Suite("W4 巴尔末系")
struct BalmerModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }
    private func rH(_ cs: ConstantsSet) throws -> Double {
        try BalmerMath.rydbergH(rInf: cs.value("R_inf"),
                                me: cs.value("m_e"), mp: cs.value("m_p"))
    }

    @Test("check1：R_H 与 CODATA 10967758.340 m⁻¹ 一致（rel < 1e-5）")
    func rydbergHydrogen() throws {
        let rh = try rH(constants())
        #expect(abs(rh - 10967758.340) / 10967758.340 < 1e-5, "R_H = \(rh)")
    }

    @Test("check2：巴尔末四线与真空观测标称值 rel < 1e-3")
    func balmerLines() throws {
        let rh = try rH(constants())
        let obs: [(n2: Int, lam: Double)] = [
            (3, 656.469), (4, 486.273), (5, 434.173), (6, 410.293),
        ]
        for (n2, o) in obs {
            let lam = BalmerMath.wavelength(n1: 2, n2: n2, rH: rh) * 1e9
            let rel = abs(lam - o) / o
            #expect(rel < 1e-3, "n2=\(n2): \(lam) nm vs \(o) nm, rel=\(rel)")
        }
    }

    @Test("能级：E_n = −13.598434/n² eV（n=1 基态 −13.598434）")
    func energyLevels() {
        #expect(abs(BalmerMath.energyLevel(1) + 13.598434) < 1e-12)
        #expect(abs(BalmerMath.energyLevel(2) + 13.598434 / 4) < 1e-12)
        #expect(abs(BalmerMath.energyLevel(6) + 13.598434 / 36) < 1e-12)
    }

    @Test("compute 输出结构：谱线散点双系列 + 能级图 + 摘要 + 理论卡")
    func computeStructure() async throws {
        let module = BalmerModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        guard case .scatter(let s) = result.charts[0],
              case .levelDiagram(let lv) = result.charts[1] else {
            Issue.record("应为 scatter + levelDiagram"); return
        }
        #expect(s.series.count == 2)
        #expect(s.series[0].points.count == 4)   // 4 条谱线计算值
        #expect(s.series[1].points.count == 4)   // 4 个观测值
        #expect(lv.levels.count == 6)
        #expect(lv.transitions.count == 3)
        #expect(result.summary.count == 5)       // R_H + 4 线
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("多值对比组：只选 Hα 时仅 1 条谱线")
    func selectionFilter() async throws {
        let module = BalmerModule()
        var values = ParamValues.defaults(for: module.params)
        values.multiCompares["lines"] = ["ha"]
        let result = try await module.compute(values, constants: try constants())
        guard case .scatter(let s) = result.charts[0] else { return }
        #expect(s.series[0].points.count == 1)
        #expect(result.summary.count == 2)       // R_H + Hα
    }

    @Test("实时档预算：compute 中位 < 16 ms")
    func computeBudget() async throws {
        let module = BalmerModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }
}
