import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W5-3 多电子原子：Slater Z_eff 全表（脚本打印值）+ NIST 电离能周期性
/// + fixtures I₁(Z) 曲线对拍。
@Suite("W5 多电子原子：Zeff 与电离能")
struct MultiElectronAtomsModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    @Test("Slater Z_eff 全表 = 脚本打印值（H..K，19 元素）")
    func slaterTable() {
        let expected: [Double] = [
            1.00, 1.70, 1.30, 1.95, 2.60, 3.25, 3.90, 4.55, 5.20, 5.85,
            2.20, 2.85, 3.50, 4.15, 4.80, 5.45, 6.10, 6.75, 2.20,
        ]
        let actual = MultiElectronAtomsMath.zeffList
        #expect(actual.count == 19)
        for i in actual.indices {
            #expect(abs(actual[i] - expected[i]) < 0.005,
                    "\(MultiElectronAtomsMath.elements[i].symbol) Z_eff = \(actual[i])")
        }
    }

    @Test("壳层锯齿：碱金属开新壳 Z_eff 骤降（Ne 5.85 → Na 2.20 → K 2.20）")
    func shellDrops() {
        let z = MultiElectronAtomsMath.zeffList
        #expect(z[10] < z[9] - 3, "Na 应骤降")
        #expect(z[18] < z[17] - 3, "K 应骤降")
        #expect(abs(z[18] - z[10]) < 1e-9, "Na 与 K 同为 2.20（Slater 合并组口径）")
    }

    @Test("电离能周期性：稀有气体峰 / 碱金属谷（NIST 打印值）")
    func ionizationPeriodicity() {
        let i1 = MultiElectronAtomsMath.ionizationEnergies
        #expect(i1.count == 20)
        #expect(abs(i1[1] - 24.59) < 0.01 && abs(i1[9] - 21.57) < 0.01
                && abs(i1[17] - 15.76) < 0.01, "He/Ne/Ar 峰")
        #expect(abs(i1[2] - 5.39) < 0.01 && abs(i1[10] - 5.14) < 0.01
                && abs(i1[18] - 4.34) < 0.01, "Li/Na/K 谷")
        // 每周期内峰谷结构：He > Li、Ne > Na、Ar > K
        #expect(i1[1] > i1[2] && i1[9] > i1[10] && i1[17] > i1[18])
    }

    @Test("fixture 曲线对拍：I₁(Z) 20 点 + 峰/谷标记系列")
    func curvesMatchFixtures() throws {
        let fx = try ModuleFixture.load("多电子原子_Zeff与电离能__v2022")
        let (x, y) = try fx.line(1, 0, index: 0)
        #expect(x.count == 20)
        for i in x.indices { #expect(x[i] == Double(i + 1), "Z = 1..20") }
        expectPointwiseClose(y, MultiElectronAtomsMath.ionizationEnergies,
                             tolerance: 1e-6, "I₁(Z)")
        // 峰（He/Ne/Ar）与谷（Li/Na/K）系列
        let (xn, yn) = try fx.line(1, 0, index: 1)
        #expect(xn == [2, 10, 18])
        #expect(abs(yn[0] - 24.587) < 1e-3 && abs(yn[1] - 21.565) < 1e-3
                && abs(yn[2] - 15.760) < 1e-3)
        let (xa, _) = try fx.line(1, 0, index: 2)
        #expect(xa == [3, 11, 19])
    }

    @Test("compute 输出结构：2 图（19 点 Z_eff + 20+3+3 点 I₁）+ 摘要 + 理论卡")
    func computeStructure() async throws {
        let module = MultiElectronAtomsModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        guard case .lineSeries(let c0) = result.charts[0],
              case .lineSeries(let c1) = result.charts[1] else {
            Issue.record("应为 lineSeries"); return
        }
        #expect(c0.series[0].points.count == 19)
        #expect(c1.series.count == 3)
        #expect(c1.series[0].points.count == 20)
        #expect(c1.series[1].points.count == 3)
        #expect(c1.series[2].points.count == 3)
        #expect(c0.referenceLines.count == 1)   // 选中元素标记
        #expect(result.summary.count == 3)
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("实时档预算：compute 中位 < 16 ms")
    func computeBudget() async throws {
        let module = MultiElectronAtomsModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }
}
