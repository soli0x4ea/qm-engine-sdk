import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W6 笔记 39《反粒子与狄拉克海》：对产生阈值（含核反冲）与狄拉克海洞。
/// fixtures：反粒子与狄拉克海_对产生阈值__v2022.json（阈值带 1.0219979 MeV）、
/// 反粒子与狄拉克海_狄拉克海洞__v2022.json（正负能带各 1.3 MeV）。
@Suite("W6 反粒子与狄拉克海")
struct DiracSeaModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    // MARK: - 物理律（对产生阈值）

    @Test("物理律：无反冲阈值 E_th = 2 m_e c² ≈ 1.022 MeV（与常数精确一致）")
    func thresholdLaw() throws {
        let c = try constants()
        let eTh = try DiracSeaMath.pairThreshold(constants: c)
        let meC2 = try c.value("me_c2_MeV")
        #expect(abs(eTh - 2 * meC2) < 1e-12, "E_th = 2mc²")
        #expect(abs(eTh - 1.02199790) < 1e-8, "E_th ≈ 1.022 MeV（已知物理值）")
    }

    @Test("物理律：核反冲抬高阈值且修正 ∝ 1/A——恒正、随 A 单调下降、重核极限回到 2mc²")
    func recoilLaw() throws {
        let c = try constants()
        let me = try c.value("m_e")
        let mp = try c.value("m_p")
        let eTh = try DiracSeaMath.pairThreshold(constants: c)

        // A = 1：E_min = 2mc²(1 + m_e/M_N)，修正 = 2mc²·(m_e/m_p) ≈ 5.6×10⁻⁴ MeV
        // （阈值由 m_e·c²/J→MeV 链路换算，双精度下相对残差 ~3×10⁻¹²——为舍入级）
        let recoil1 = try DiracSeaMath.pairThresholdRecoil(constants: c, nucleusMass: mp)
        let expectedRecoil = eTh * me / mp
        #expect(recoil1 > eTh, "动量守恒 ⇒ 阈值严格高于 2mc²")
        #expect(abs(recoil1 - eTh - expectedRecoil) < 5e-12 * eTh,
                "反冲修正 = 2mc²·(m_e/M_N)（三体阈值的解析结构）")
        // ∝ 1/A：A = 10 的修正恰为 A = 1 的 1/10（差值含灾难性消去，相对精度 ~1e-8）
        let recoil10 = try DiracSeaMath.pairThresholdRecoil(constants: c, nucleusMass: 10 * mp)
        #expect(abs((recoil10 - eTh) / expectedRecoil - 0.1) < 1e-7)
        // A = 250 重核：修正缩小 250 倍，趋回 2mc²
        let recoil250 = try DiracSeaMath.pairThresholdRecoil(constants: c, nucleusMass: 250 * mp)
        #expect(recoil250 < recoil10 && recoil250 < eTh + 1e-5,
                "A→∞ 极限回到无反冲阈值")
        // 守恒律判定与符号：修正恒正
        #expect((recoil1 - eTh) > 0 && (recoil250 - eTh) > 0)
    }

    // MARK: - 物理律（狄拉克海洞）

    @Test("物理律：海面占据态全部为负能（y < 0）、落入示意带内；固定种子 N=90 确定性")
    func seaConfigurationLaw() throws {
        let eBound = 1.30
        let sea = DiracSeaMath.seaPoints(count: 90, eBound: eBound)
        #expect(sea.count == 90)
        for p in sea {
            #expect(p.y < 0, "狄拉克海：负能态（E = \(p.y) MeV < 0）")
            #expect(p.y > -eBound && p.y < -0.05, "落于负能连续谱示意区")
            #expect(p.x > -0.9 && p.x < 0.9)
        }
        // 固定种子确定性：同种子同构型（SplitMix64，规范 docs/RNG_SEED_POLICY.md）
        let again = DiracSeaMath.seaPoints(count: 90, eBound: eBound)
        #expect(zip(sea, again).allSatisfy { $0.x == $1.x && $0.y == $1.y })
        // 不同种子给出不同构型（统计意义上独立）
        let other = DiracSeaMath.seaPoints(count: 90, seed: 42, eBound: eBound)
        #expect(zip(sea, other).contains { $0.x != $1.x || $0.y != $1.y })
        // 点数可变：N = 30
        #expect(DiracSeaMath.seaPoints(count: 30, eBound: eBound).count == 30)
    }

    @Test("物理律：洞的量子数翻转——缺失 E = −mc² 电子 ⇒ 洞能级 |E₀| = mc²、电荷 +e")
    func holeQuantumNumbersLaw() async throws {
        let c = try constants()
        let e0 = try DiracSeaMath.restEnergy(constants: c)
        #expect(abs(e0 - 0.51099895) < 1e-7, "|E₀| = m_e c² ≈ 0.511 MeV")
        // 模块输出：洞标记在 y = −mc²，正电子标记在 y = +mc²（能量镜像）
        let module = DiracSeaHoleModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: c)
        guard case .schematic(let d) = result.charts[0] else {
            Issue.record("应为 schematic"); return
        }
        let hole = d.markers.first { $0.label?.contains("洞") == true }
        let positron = d.markers.first { $0.label?.contains("正电子") == true }
        #expect(abs((hole?.y ?? 0) + e0) < 1e-12, "洞标记于 y = −mc²")
        #expect(abs((positron?.y ?? 0) - e0) < 1e-12, "正电子标记于 y = +mc²（符号翻转）")
        #expect(abs((hole?.y ?? 0) + (positron?.y ?? 1)) < 1e-12,
                "洞 ↔ 正电子能量镜像：|−E₀| = |+E₀|")
    }

    // MARK: - fixtures 对拍

    @Test("fixtures：阈值带高度 = 2mc²；海图正负能带各 1.3 MeV")
    func bandsMatchFixtures() throws {
        let c = try constants()
        let fxThreshold = try ModuleFixture.load("反粒子与狄拉克海_对产生阈值__v2022")
        let bars = try fxThreshold.bars(0, 0)
        #expect(bars.count == 1)
        let threshold = try DiracSeaMath.pairThreshold(constants: c)
        #expect(abs(bars[0].height - threshold) < 1e-8,
                "阈值带高度 \(bars[0].height) = 2mc²")

        let fxSea = try ModuleFixture.load("反粒子与狄拉克海_狄拉克海洞__v2022")
        let bands = try fxSea.bars(0, 0)
        #expect(bands.count == 2, "负能海带 + 正能连续谱带")
        for b in bands { #expect(abs(b.height - 1.3) < 1e-8) }
    }

    // MARK: - compute 结构与预算

    @Test("compute 结构与预算：两模块实时档（阈值示意带 + 海面 90 点）")
    func computeStructureAndBudget() async throws {
        let c = try constants()

        let threshold = PairProductionModule()
        let r1 = try await threshold.compute(
            ParamValues.defaults(for: threshold.params), constants: c)
        #expect(r1.charts.count == 1)
        guard case .schematic(let s1) = r1.charts[0] else {
            Issue.record("阈值模块应为 schematic"); return
        }
        #expect(s1.bands.count == 1 && s1.callouts.count == 3)
        #expect(r1.summary.count == 4)
        #expect(r1.theory?.formulas.count == 4)
        try await expectComputeUnderBudget(
            module: threshold, values: ParamValues.defaults(for: threshold.params),
            constants: c, budgetMillis: 16)

        let sea = DiracSeaHoleModule()
        let r2 = try await sea.compute(
            ParamValues.defaults(for: sea.params), constants: c)
        #expect(r2.charts.count == 1)
        guard case .schematic(let s2) = r2.charts[0] else {
            Issue.record("海模块应为 schematic"); return
        }
        #expect(s2.bands.count == 2)
        #expect(s2.markers.count == 93, "90 海点 + 洞 + 电子 + 正电子")
        #expect(r2.summary.count == 5)
        #expect(r2.theory?.formulas.count == 4)
        try await expectComputeUnderBudget(
            module: sea, values: ParamValues.defaults(for: sea.params),
            constants: c, budgetMillis: 16)
    }
}
