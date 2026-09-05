import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W6 笔记 33《量子隐形传态与密集编码》：Werner 信道保真度 / 密集编码容量。
/// fixtures：量子隐形传态_保真度__v2022.json（F(p) 401 点 + 两条参考线）、
/// 量子隐形传态_密集编码容量__v2022.json（I_total(p) 401 点 + Holevo/容量参考线）。
@Suite("W6 量子隐形传态与密集编码")
struct TeleportationModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    // MARK: - 物理律

    @Test("物理律：F(0)=1/2、F(1/3)=2/3（经典极限恰在可分边界）、F(1)=1；单调递增")
    func fidelityLaw() {
        #expect(abs(WernerChannelMath.fidelity(p: 0) - 0.5) < 1e-15, "最大混合信道 F = 1/2")
        #expect(abs(WernerChannelMath.fidelity(p: 1.0 / 3.0) - 2.0 / 3.0) < 1e-15,
                "可分边界 p=1/3 处 F 恰为经典极限 2/3")
        #expect(abs(WernerChannelMath.fidelity(p: 1) - 1) < 1e-15, "纯贝尔信道完美传态 F = 1")
        // 单调性：偏导 1/2 > 0
        for p in stride(from: 0.0, through: 0.99, by: 0.01) {
            #expect(WernerChannelMath.fidelity(p: p + 0.01) > WernerChannelMath.fidelity(p: p))
        }
        // 量子优势边界：F > 2/3 ⟺ p > 1/3
        #expect(WernerChannelMath.fidelity(p: 0.34) > 2.0 / 3.0)
        #expect(WernerChannelMath.fidelity(p: 0.32) < 2.0 / 3.0)
    }

    @Test("物理律：von Neumann 熵 S(0)=2 bit、S(1)=0；单态分数 f(1/3)=1/2")
    func entropyLaw() {
        #expect(abs(WernerChannelMath.vonNeumannEntropy(p: 0) - 2) < 1e-15,
                "p=0 为最大混合态 I/4，S = 2 bit")
        #expect(abs(WernerChannelMath.vonNeumannEntropy(p: 1)) < 1e-15,
                "p=1 为纯贝尔态，S = 0")
        #expect(abs(WernerChannelMath.singletFraction(p: 1.0 / 3.0) - 0.5) < 1e-15)
        // 0 ≤ S ≤ 2 且随 p 单调下降（纠缠纯化降低熵）
        for p in stride(from: 0.0, through: 0.99, by: 0.05) {
            let s = WernerChannelMath.vonNeumannEntropy(p: p)
            #expect(s >= 0 && s <= 2)
            #expect(WernerChannelMath.vonNeumannEntropy(p: p + 0.05)
                    < WernerChannelMath.vonNeumannEntropy(p: p))
        }
    }

    @Test("物理律：Holevo 界 I ≥ 1 恒成立；增益阈值 p* ≈ 0.748（S(ρ)=1，高于可分边界 1/3）；I(0)=1、I(1)=2")
    func denseCodingLaw() {
        #expect(abs(WernerChannelMath.denseCodingCapacity(p: 0) - 1) < 1e-15,
                "无纠缠资源：退化为 Holevo 界 1 bit")
        #expect(abs(WernerChannelMath.denseCodingCapacity(p: 1) - 2) < 1e-15,
                "纯贝尔资源：1 ebit + 1 qubit = 2 bit")
        for p in stride(from: 0.0, through: 1.0 / 3.0, by: 1.0 / 30.0) {
            #expect(abs(WernerChannelMath.denseCodingCapacity(p: p) - 1) < 1e-15,
                    "p=\(p)：可分资源，纠缠增益 C_extra < 0 被 max(0,·) 截断")
        }
        // 纠缠 ≠ 可用：p ∈ (1/3, p*) 资源态纠缠，但 S(ρ) > 1 ⇒ 容量仍钉在 Holevo 界
        for p in [0.4, 0.5, 0.6, 0.7] {
            #expect(abs(WernerChannelMath.denseCodingCapacity(p: p) - 1) < 1e-15,
                    "p=\(p)：S(ρ)=\(String(format: "%.4f", WernerChannelMath.vonNeumannEntropy(p: p))) > 1，无密集编码增益")
        }
        for p in [0.75, 0.85, 0.95] {
            let i = WernerChannelMath.denseCodingCapacity(p: p)
            #expect(i > 1 && i <= 2)
        }
        // 增益阈值自洽：S(p*) = 1（二分根），且 p* ∈ (1/3, 1)
        let pStar = WernerChannelMath.gainThreshold
        #expect(abs(WernerChannelMath.vonNeumannEntropy(p: pStar) - 1) < 1e-12)
        #expect(abs(pStar - 0.7476) < 1e-3, "文献值 p* ≈ 0.7476，实测 \(pStar)")
        #expect(abs(WernerChannelMath.denseCodingCapacity(p: pStar + 0.01)
                    - WernerChannelMath.denseCodingCapacity(p: pStar - 0.01)) > 0.01,
                "跨过 p* 容量脱离 Holevo 界")
    }

    // MARK: - fixtures 对拍

    @Test("fixtures：F(p) 401 点逐点对拍（< 1e-8）+ 经典极限 / 完美传态参考线")
    func fidelityMatchesFixtures() throws {
        let fx = try ModuleFixture.load("量子隐形传态_保真度__v2022")
        let (ps, fF) = try fx.line(0, 0, label: "F(p)")
        #expect(ps.count == 401)
        expectPointwiseClose(ps.map { WernerChannelMath.fidelity(p: $0) }, fF,
                             tolerance: 1e-8, "F(p) = (p+1)/2")

        let (_, classical) = try fx.line(0, 0, label: "classical limit")
        #expect(abs(classical[0] - 2.0 / 3.0) < 1e-8 && classical[0] == classical[1])
        let (_, perfect) = try fx.line(0, 0, label: "perfect")
        #expect(abs(perfect[0] - 1) < 1e-8)
    }

    @Test("fixtures：I_total(p) 401 点逐点对拍（< 1e-8）+ Holevo / 最大容量参考线")
    func denseCodingMatchesFixtures() throws {
        let fx = try ModuleFixture.load("量子隐形传态_密集编码容量__v2022")
        let (ps, iF) = try fx.line(0, 0, label: "total")
        #expect(ps.count == 401)
        expectPointwiseClose(ps.map { WernerChannelMath.denseCodingCapacity(p: $0) }, iF,
                             tolerance: 1e-8, "I_total(p) = 1 + max(0, 1 − S(ρ))")

        let (_, holevo) = try fx.line(0, 0, label: "Holevo")
        #expect(abs(holevo[0] - 1) < 1e-8)
        let (_, maxCap) = try fx.line(0, 0, label: "max capacity")
        #expect(abs(maxCap[0] - 2) < 1e-8)
    }

    // MARK: - compute 结构与预算

    @Test("compute 结构：两模块各 1 图 2 系列（401 点）+ 摘要 + 理论卡；实时档预算")
    func computeStructureAndBudget() async throws {
        let modules: [any SimModule] = [TeleportationFidelityModule(), DenseCodingModule()]
        for module in modules {
            let result = try await module.compute(
                ParamValues.defaults(for: module.params), constants: try constants())
            #expect(result.charts.count == 1)
            guard case .lineSeries(let c) = result.charts[0] else {
                Issue.record("\(module.meta.id) 应为 lineSeries"); return
            }
            #expect(c.series[0].points.count == 401)
            #expect(result.summary.count == 4)
            let expectedFormulas = module.meta.id.contains("密集编码") ? 5 : 4
            #expect(result.theory?.formulas.count == expectedFormulas)

            try await expectComputeUnderBudget(
                module: module, values: ParamValues.defaults(for: module.params),
                constants: try constants(), budgetMillis: 16)
        }
    }
}
