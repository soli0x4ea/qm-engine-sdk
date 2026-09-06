import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W6 笔记 37《克莱因-戈登方程》：色散关系 E = ±√(p²c² + m²c⁴)。
/// fixtures：克莱因-戈登方程_色散关系__v2022.json（三线各 512 点，m₀ = m_e c²）。
@Suite("W6 克莱因-戈登色散关系")
struct KleinGordonModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }
    private static let m0 = 0.51099895069   // m_e c²（MeV）

    // MARK: - 物理律

    @Test("物理律：p = 0 质量阈值 E = ±mc²；正负能隙 2mc²")
    func massThresholdLaw() {
        #expect(KGMath.energyPositive(0, m0: Self.m0) == Self.m0, "E₊(0) = mc²（浮点精确）")
        #expect(KGMath.energyNegative(0, m0: Self.m0) == -Self.m0)
        // 正负能支在 p=0 处的间隙 = 2mc²（狄拉克海的伏笔）
        #expect(abs(KGMath.energyPositive(0, m0: Self.m0)
                    - KGMath.energyNegative(0, m0: Self.m0) - 2 * Self.m0) < 1e-15)
    }

    @Test("物理律：正能支随 |p| 单调抬升；E₋ = −E₊ 处处成立")
    func symmetryAndMonotonicityLaw() {
        let pcs = stride(from: -2.0, through: 2.0, by: 0.17).map { $0 }
        for i in pcs.indices.dropFirst() {
            if abs(pcs[i]) > abs(pcs[i - 1]) {
                #expect(KGMath.energyPositive(pcs[i], m0: Self.m0)
                        > KGMath.energyPositive(pcs[i - 1], m0: Self.m0))
            }
            #expect(abs(KGMath.energyPositive(pcs[i], m0: Self.m0)
                        + KGMath.energyNegative(pcs[i], m0: Self.m0)) < 1e-15,
                    "E₋(p) = −E₊(p)")
        }
    }

    @Test("物理律：非相对论极限——E₊ → mc² + p²/2m（误差 ∝ (pc/mc²)⁴）；高动量 NR 高估")
    func nonrelativisticLimitLaw() {
        // pc = mc²/100 时：E₊ − NR ≈ −mc²·(pc/mc²)⁴/8 = −1.25×10⁻⁹·mc²
        let small = Self.m0 / 100
        #expect(abs(KGMath.energyPositive(small, m0: Self.m0)
                    - KGMath.energyNR(small, m0: Self.m0)) < 1.5e-9 * Self.m0,
                "小动量下 NR 近似贴合到 (pc/mc²)⁴ 阶")
        // 误差随动量四次方增长：pc × 10 → 误差 × 10⁴
        let e1 = abs(KGMath.energyPositive(small, m0: Self.m0)
                     - KGMath.energyNR(small, m0: Self.m0))
        let e2 = abs(KGMath.energyPositive(10 * small, m0: Self.m0)
                     - KGMath.energyNR(10 * small, m0: Self.m0))
        #expect(e2 / e1 > 5_000 && e2 / e1 < 15_000, "误差比 \(e2 / e1) ≈ 10⁴")
        // NR 抛物线在高动量端高估能量（无光速上限）
        let hi = 10 * Self.m0
        #expect(KGMath.energyNR(hi, m0: Self.m0) > KGMath.energyPositive(hi, m0: Self.m0))
        // 极端相对论质量占比 (mc²/E)²：pc = 10mc² → 1/101
        #expect(abs(KGMath.massFraction(hi, m0: Self.m0) - 1.0 / 101.0) < 1e-15)
    }

    // MARK: - fixtures 对拍

    @Test("fixtures：正能 / 负能 / 非相对论三线 512 点逐点（< 1e-8）")
    func matchesFixtures() throws {
        let fx = try ModuleFixture.load("克莱因-戈登方程_色散关系__v2022")
        let (xPos, yPos) = try fx.line(0, 0, label: "E = +")
        let (xNeg, yNeg) = try fx.line(0, 0, label: "E = -")
        let (xNR, yNR) = try fx.line(0, 0, label: "NR")
        #expect(xPos.count == 512 && xNeg.count == 512 && xNR.count == 512)
        #expect(xPos == xNR, "三线共用同一 p·c 网格")

        expectPointwiseClose(xPos.map { KGMath.energyPositive($0, m0: Self.m0) },
                              yPos, tolerance: 1e-8, "E₊(p)")
        expectPointwiseClose(xNeg.map { KGMath.energyNegative($0, m0: Self.m0) },
                              yNeg, tolerance: 1e-8, "E₋(p)")
        expectPointwiseClose(xNR.map { KGMath.energyNR($0, m0: Self.m0) },
                              yNR, tolerance: 1e-8, "NR 近似")
    }

    // MARK: - compute 结构与预算

    @Test("compute 结构：1 图 3 系列（500 点）+ ±mc² 参考线 + 摘要 + 理论卡；实时档预算")
    func computeStructureAndBudget() async throws {
        let module = KleinGordonDispersionModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 1)
        guard case .lineSeries(let c) = result.charts[0] else {
            Issue.record("应为 lineSeries"); return
        }
        #expect(c.series.count == 3)
        // W13h #37：NR 抛物线出画截断（y ≤ 1.02·maxE₊）——E₊/E₋ 仍 500 点
        #expect(c.series[0].points.count == 500)
        #expect(c.series[1].points.count == 500)
        let nr = c.series[2].points
        #expect(nr.count < 500 && nr.count > 100, "NR 截断后 \(nr.count) 点")
        let topY = c.series[0].points.map(\.y).max()! * 1.02
        #expect(nr.allSatisfy { $0.y <= topY }, "NR 全部在画窗内")
        #expect(c.referenceLines.count == 3)
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 4)

        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }
}
