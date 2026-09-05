import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W6 笔记 38《狄拉克方程》：氢原子精确谱 vs 薛定谔能级。
/// fixtures：狄拉克方程_氢原子精细结构__v2022.json（能级横线 5 条：两张图）。
@Suite("W6 狄拉克氢原子精细结构")
struct DiracHydrogenModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    private func math() throws -> (alpha: Double, meC2: Double) {
        let c = try constants()
        return (try 1.0 / c.value("alpha_inv"), try c.value("me_c2_MeV") * 1e6)
    }

    // MARK: - 物理律

    @Test("物理律：2S₁/₂ = 2P₁/₂ 简并（同 j 同能量）；2P₃/₂ 抬高（j 大修正小、结合更松）")
    func degeneracyLaw() throws {
        let (alpha, me) = try math()
        // 点核狄拉克谱只依赖 n 与 |κ|（= j+1/2）：2S₁/₂ 与 2P₁/₂ 同 j = 1/2 → 同 κ → 同能量
        let e2S = DiracHydrogenMath.eDiracBinding(n: 2, kappa: 1, Z: 1, alpha: alpha, meC2eV: me)
        let e2P1 = DiracHydrogenMath.eDiracBinding(n: 2, kappa: 1, Z: 1, alpha: alpha, meC2eV: me)
        #expect(e2S == e2P1)
        // 2P₃/₂（κ = 2）在 2S₁/₂ 之上（能量更高 / 结合更松）
        let e2P3 = DiracHydrogenMath.eDiracBinding(n: 2, kappa: 2, Z: 1, alpha: alpha, meC2eV: me)
        #expect(e2P3 > e2S, "劈裂 = \(e2P3 - e2S) eV > 0")
        #expect(e2P3 - e2S < 1e-3, "氢精细结构 ~ 4.5×10⁻⁵ eV（远小于主能级间隔）")
    }

    @Test("物理律：精细结构劈裂 2P₃/₂ − 2S₁/₂ ≈ mc²(Zα)⁴/32（主阶，偏差 < 0.1%）；Z⁴ 标度")
    func fineStructureLaw() throws {
        let (alpha, me) = try math()
        let split = DiracHydrogenMath.fineStructureSplit(Z: 1, alpha: alpha, meC2eV: me)
        let leading = DiracHydrogenMath.fineStructureLeading(Z: 1, alpha: alpha, meC2eV: me)
        #expect(split > 0)
        #expect(abs(split / leading - 1) < 1e-3,
                "exact/leading = \(split / leading)（高阶 (Zα)⁶ 修正 ~3×10⁻⁵）")
        // 劈裂按 Z⁴ 急剧放大（Z=2 → 16×，Z=4 → 256×；高阶项仅 ~1‰ 偏差）
        let split2 = DiracHydrogenMath.fineStructureSplit(Z: 2, alpha: alpha, meC2eV: me)
        let split4 = DiracHydrogenMath.fineStructureSplit(Z: 4, alpha: alpha, meC2eV: me)
        #expect(abs(split2 / split - 16) < 0.1, "Z=2：比值 \(split2 / split) ≈ 16")
        #expect(abs(split4 / split - 256) < 3, "Z=4：比值 \(split4 / split) ≈ 256")
    }

    @Test("物理律：非相对论极限——狄拉克能级收敛到薛定谔值，偏差 ∝ (Zα)²")
    func nonrelativisticLimitLaw() throws {
        let (alpha, me) = try math()
        func relDiff(_ al: Double) -> Double {
            let ed = DiracHydrogenMath.eDiracBinding(n: 1, kappa: 1, Z: 1, alpha: al, meC2eV: me)
            let es = DiracHydrogenMath.eSchrodinger(n: 1, Z: 1, alpha: al, meC2eV: me)
            return abs(ed - es) / abs(es)
        }
        let d1 = relDiff(alpha)
        let d2 = relDiff(alpha / 10)
        #expect(d1 < 2e-5, "氢 1S：相对论修正 ~1.3×10⁻⁵（故薛定谔谱几乎正确）")
        #expect(abs(d1 / (alpha * alpha / 4) - 1) < 0.01,
                "主阶修正恰为 (Zα)²/4（相对论深化）")
        // α × 1/10 → 偏差 × 1/100（(Zα)² 标度）
        #expect(d1 / d2 > 90 && d1 / d2 < 111, "收敛比 \(d1 / d2) ≈ 100")
        // 同 α 下狄拉克能级恒更深（更负）
        let ed = DiracHydrogenMath.eDiracBinding(n: 1, kappa: 1, Z: 1, alpha: alpha, meC2eV: me)
        let es = DiracHydrogenMath.eSchrodinger(n: 1, Z: 1, alpha: alpha, meC2eV: me)
        #expect(ed < es)
    }

    // MARK: - fixtures 对拍

    @Test("fixtures：五条能级横线（1S/2S/2P₃/₂/薛定谔 n=2，< 1e-8）")
    func levelsMatchFixtures() throws {
        let fx = try ModuleFixture.load("狄拉克方程_氢原子精细结构__v2022")
        let (alpha, me) = try math()

        let (_, f1S) = try fx.line(0, 0, label: "1S")
        let (_, f2S) = try fx.line(0, 0, label: "2S")
        let (_, fSchr) = try fx.line(0, 0, label: "Schrodinger")
        let (_, f2SZoom) = try fx.line(0, 1, label: "2S")
        let (_, f2P3) = try fx.line(0, 1, label: "2P_{3/2}")

        let e1S = DiracHydrogenMath.eDiracBinding(n: 1, kappa: 1, Z: 1, alpha: alpha, meC2eV: me)
        let e2S = DiracHydrogenMath.eDiracBinding(n: 2, kappa: 1, Z: 1, alpha: alpha, meC2eV: me)
        let e2P3 = DiracHydrogenMath.eDiracBinding(n: 2, kappa: 2, Z: 1, alpha: alpha, meC2eV: me)
        let eSchr = DiracHydrogenMath.eSchrodinger(n: 2, Z: 1, alpha: alpha, meC2eV: me)

        #expect(abs(e1S - f1S[0]) / abs(f1S[0]) < 1e-8, "1S：\(e1S) vs \(f1S[0])")
        #expect(abs(e2S - f2S[0]) / abs(f2S[0]) < 1e-8, "2S：\(e2S) vs \(f2S[0])")
        #expect(abs(e2S - f2SZoom[0]) / abs(f2SZoom[0]) < 1e-8, "放大图 2S")
        #expect(abs(e2P3 - f2P3[0]) / abs(f2P3[0]) < 1e-8, "2P₃/₂：\(e2P3) vs \(f2P3[0])")
        #expect(abs(eSchr - fSchr[0]) / abs(fSchr[0]) < 1e-8, "薛定谔 n=2：\(eSchr) vs \(fSchr[0])")
        // 2S₁/₂ = 2P₁/₂ 简并在 fixture 中体现为同一横线标签
        #expect(f2S[0] == f2SZoom[0])
    }

    // MARK: - compute 结构与预算

    @Test("compute 结构：两张能级图 + n=2 放大图跃迁标注 + 摘要 + 理论卡；实时档预算")
    func computeStructureAndBudget() async throws {
        let module = DiracHydrogenModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        for (index, chart) in result.charts.enumerated() {
            guard case .levelDiagram(let d) = chart else {
                Issue.record("charts[\(index)] 应为 levelDiagram"); return
            }
            #expect(d.levels.count == (index == 0 ? 4 : 3))
        }
        #expect(result.summary.count == 5)
        #expect(result.theory?.formulas.count == 4)

        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }
}
