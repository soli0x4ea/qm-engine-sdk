import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W5-4 分子结构：Morse 势端点/极值（脚本打印值）+ Hückel 能级
/// + fixtures Morse 曲线对拍。
@Suite("W5 分子结构：H₂ 与苯")
struct MolecularStructureModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    /// 脚本默认参数：D_e = 4.52 eV，a = 1.94 Å⁻¹，R_e = 0.741 Å。
    private var defaults: (de: Double, a: Double, re: Double) { (4.52, 1.94, 0.741) }

    @Test("Morse 势：E(R_e) = −D_e，E(4 Å) ≈ −0.016 eV（脚本打印值）")
    func morseEndpoints() {
        let p = defaults
        #expect(abs(MolecularStructureMath.morse(p.re, de: p.de, a: p.a, re: p.re) + 4.52) < 1e-9)
        #expect(abs(MolecularStructureMath.morse(4.0, de: p.de, a: p.a, re: p.re) + 0.016) < 1e-3)
        // 井口收敛：R→∞ 上界 0
        #expect(MolecularStructureMath.morse(20, de: p.de, a: p.a, re: p.re) < 0)
        #expect(MolecularStructureMath.morse(20, de: p.de, a: p.a, re: p.re) > -1e-12)
    }

    @Test("Morse 井：R_e 为唯一极小、两侧抬升")
    func morseWell() {
        let p = defaults
        let vRe = MolecularStructureMath.morse(p.re, de: p.de, a: p.a, re: p.re)
        for d in [0.05, 0.15, 0.3, 0.5] {
            let vl = MolecularStructureMath.morse(p.re - d, de: p.de, a: p.a, re: p.re)
            let vr = MolecularStructureMath.morse(p.re + d, de: p.de, a: p.a, re: p.re)
            #expect(vl > vRe && vr > vRe)
        }
    }

    @Test("苯 Hückel 能级：α+2βcos(2πk/6)，HOMO/LUMO/隙（脚本打印值）")
    func huckel() {
        let beta = -2.5
        let eps = MolecularStructureMath.huckelLevels(beta: beta)
        #expect(eps.count == 6)
        let expected: [Double] = [-5.0, -2.5, 2.5, 5.0, 2.5, -2.5]
        for i in eps.indices { #expect(abs(eps[i] - expected[i]) < 1e-12, "k=\(i)") }
        // 二重简并：k=1,5 与 k=2,4
        #expect(abs(eps[1] - eps[5]) < 1e-12 && abs(eps[2] - eps[4]) < 1e-12)
        let gap = MolecularStructureMath.homoLumoGap(beta: beta)
        #expect(abs(gap.homo + 2.5) < 1e-12, "HOMO = α+β")
        #expect(abs(gap.lumo - 2.5) < 1e-12, "LUMO = α−β")
        #expect(abs(gap.gap - 5.0) < 1e-12, "π-π* 隙 = 2|β|")
    }

    @Test("fixture 曲线对拍：Morse V(R)（400 点全量）")
    func curvesMatchFixtures() throws {
        let fx = try ModuleFixture.load("分子结构_H2与苯__v2022")
        let (x, y) = try fx.line(0, 0, index: 0)
        #expect(x.count == 400)
        let p = defaults
        // Morse 含 e^(−2aR) 指数，x 舍入传导略大 → 容差 1e-7
        expectPointwiseClose(x.map { MolecularStructureMath.morse($0, de: p.de, a: p.a, re: p.re) },
                             y, tolerance: 1e-7, "Morse V(R)")
    }

    @Test("compute 输出结构：Morse 400 点 + 能级图 6 级 + 摘要 + 理论卡")
    func computeStructure() async throws {
        let module = MolecularStructureModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        guard case .lineSeries(let c0) = result.charts[0],
              case .levelDiagram(let c1) = result.charts[1] else {
            Issue.record("应为 lineSeries + levelDiagram"); return
        }
        #expect(c0.series[0].points.count == 400)
        #expect(c0.referenceLines.count == 2)   // R_e 与 −D_e
        #expect(c1.levels.count == 6)
        // 占据：k=0,1,5（6 π 电子）；能级升序后校验
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("实时档预算：compute 中位 < 16 ms")
    func computeBudget() async throws {
        let module = MolecularStructureModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }
}
