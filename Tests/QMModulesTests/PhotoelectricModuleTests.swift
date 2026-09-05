import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W3-3 光电效应模块：fixtures 对拍（光电效应_截止电压线性__v2022.json）。
///
/// fixture x 存 9 位有效数字，在 V₀→0 端（y ~ 5×10⁻³）经斜率 h/e 放大为 ~4×10⁻⁷
/// 相对误差，对拍容差取 1×10⁻⁵。
@Suite("W3 光电效应模块")
struct PhotoelectricModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    @Test("四金属截止电压曲线逐点对拍 fixtures")
    func linesMatchFixtures() throws {
        let fx = try ModuleFixture.load("光电效应_截止电压线性__v2022")
        let cs = try constants()
        let hE = PhotoelectricMath.hOverE(h: try cs.value("h"), e: try cs.value("e"))

        // 复刻脚本断言：斜率普适 h/e = 4.135667e-15 V·s
        #expect(abs(hE - 4.135667e-15) < 1e-20, "h/e = \(hE)")

        for m in PhotoelectricModule.metals {
            let (xs, ys) = try fx.line(0, 0, label: "\(m.id) (")
            expectPointwiseClose(
                xs.map { PhotoelectricMath.stoppingVoltage($0 * 1e15, W: m.work, hOverE: hE) },
                ys, tolerance: 1e-5, "\(m.id) 截止电压直线")
        }
    }

    @Test("compute 输出结构：四系列点数 / 参考线 / λ₀ 摘要")
    func computeStructure() async throws {
        let module = PhotoelectricModule()
        let values = ParamValues.defaults(for: module.params)
        let result = try await module.compute(values, constants: try constants())

        guard case .lineSeries(let line) = result.charts[0] else {
            Issue.record("charts[0] 应为 lineSeries"); return
        }
        #expect(line.series.map(\.name) == ["Na", "K", "Cs", "W"])
        // linspace(0, 1.6e15, 400) 上 V₀≥0 的点数（网格值距 ν₀ 最小 |V₀| ≈ 5e-3，
        // 无边界翻转风险，与 Python mask 逐点一致）
        let expectedCounts = ["Na": 262, "K": 261, "Cs": 273, "W": 128]
        for s in line.series {
            #expect(s.points.count == expectedCounts[s.name],
                    "\(s.name) 点数 \(s.points.count) vs \(expectedCounts[s.name] ?? -1)")
            // 物理区全部 V₀ ≥ 0
            #expect(s.points.allSatisfy { $0.y >= 0 })
        }
        // 四条 ν₀ 参考线
        #expect(line.referenceLines.count == 4)
        let cs = try constants()
        let hE = PhotoelectricMath.hOverE(h: try cs.value("h"), e: try cs.value("e"))
        let nu0Na = PhotoelectricMath.cutoffs(W: 2.28, hOverE: hE, hcEVNm: 1239.84).nu0
        let refNa = line.referenceLines.first { $0.id == "nu0_Na" }
        #expect(refNa.map { abs($0.value - nu0Na / 1e15) < 0.01 } == true)

        // 摘要：斜率 + hc + 四金属 λ₀（与 fixture 图例 λ₀=544/539/590/276 nm 一致）
        #expect(result.summary.count == 6)
        let lam0 = Dictionary(uniqueKeysWithValues:
            result.summary.filter { $0.id.hasPrefix("lam0_") }.map { (String($0.id.dropFirst(5)), $0.value) })
        #expect(lam0["Na"]?.hasPrefix("544") == true)
        #expect(lam0["K"]?.hasPrefix("539") == true)
        #expect(lam0["Cs"]?.hasPrefix("590") == true)
        #expect(lam0["W"]?.hasPrefix("276") == true)
        #expect(result.theory?.formulas.count == 3)
    }

    @Test("实时档预算：compute 中位 < 16 ms")
    func computeBudget() async throws {
        let module = PhotoelectricModule()
        try await expectComputeUnderBudget(
            module: module,
            values: ParamValues.defaults(for: module.params),
            constants: try constants(),
            budgetMillis: 16)
    }
}
