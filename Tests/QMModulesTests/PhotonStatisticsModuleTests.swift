import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W3-4 光子统计模块：fixtures 对拍（光场量子化_相干态与光子统计__v2022.json，μ = 5）。
@Suite("W3 光子统计模块")
struct PhotonStatisticsModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    @Test("泊松 / 玻色-爱因斯坦 26 柱逐柱对拍 fixtures")
    func barsMatchFixtures() throws {
        let fx = try ModuleFixture.load("光场量子化_相干态与光子统计__v2022")
        let mu = 5.0

        let coh = try fx.bars(0, 0)
        let th = try fx.bars(0, 1)
        #expect(coh.count == 26 && th.count == 26, "n_show = 25 → 各 26 柱")
        #expect(coh.map(\.x) == (0...25).map { Double($0) })

        expectPointwiseClose(
            (0...25).map { PhotonStatsMath.coherent($0, mu: mu) },
            coh.map(\.height), tolerance: 1e-8, "相干态 P(n)")
        expectPointwiseClose(
            (0...25).map { PhotonStatsMath.thermal($0, mu: mu) },
            th.map(\.height), tolerance: 1e-8, "热态 P(n)")
    }

    @Test("统计矩：归一化 / 均值 / 方差 / Mandel Q（复刻脚本打印值）")
    func momentsAndMandelQ() {
        let mu = 5.0
        // 求和域：Python 固定 nmax = 120，μ = 5 时尾部已收敛
        #expect(PhotonStatsMath.nCompute(forMu: mu) == 120)

        let nMax = PhotonStatsMath.nCompute(forMu: mu)
        let coh = PhotonStatsMath.distribution(PhotonStatsMath.coherent, mu: mu, n: nMax)
        let th = PhotonStatsMath.distribution(PhotonStatsMath.thermal, mu: mu, n: nMax)
        let (mCoh, m2Coh) = PhotonStatsMath.moments(coh)
        let (mTh, m2Th) = PhotonStatsMath.moments(th)

        #expect(abs(coh.reduce(0, +) - 1) < 1e-10, "相干态 ΣP")
        #expect(abs(th.reduce(0, +) - 1) < 1e-9, "热态 ΣP（尾部 (5/6)¹²¹ ≈ 3×10⁻¹⁰）")
        // 泊松：⟨n⟩ = Δn² = μ = 5（截断尾部可忽略）
        #expect(abs(mCoh - 5) < 1e-9)
        #expect(abs(m2Coh - mCoh * mCoh - 5) < 1e-8)
        // 热态：⟨n⟩ = 5，Δn² = μ(1+μ) = 30。n_max = 120 截断尾部 (5/6)¹²¹ ≈ 2.7e-10，
        // 均值/方差损失 ~3e-8 / ~4e-6（Python 同参数同样偏差，打印 6 位小数仍为 5.000000）
        #expect(abs(mTh - 5) < 1e-7)
        #expect(abs(m2Th - mTh * mTh - 30) < 1e-5)
        // Mandel Q（热态超泊松度）= μ = 5（含截断效应 ~7e-7）
        let q = (m2Th - mTh * mTh) / mTh - 1
        #expect(abs(q - 5) < 1e-5, "Mandel Q = \(q)")
    }

    @Test("compute 输出结构：双图 26 柱 / 摘要五卡")
    func computeStructure() async throws {
        let module = PhotonStatisticsModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())

        #expect(result.charts.count == 2)
        for (index, chart) in result.charts.enumerated() {
            guard case .bars(let bars) = chart else {
                Issue.record("charts[\(index)] 应为 bars"); return
            }
            #expect(bars.bars.count == 26, "μ=5、nshow=25 → \(index == 0 ? "相干态" : "热态") 26 柱")
        }
        #expect(result.summary.count == 5)
        #expect(result.theory?.formulas.count == 3)
    }

    @Test("高端参数 μ=50：求和域自适应扩展且仍满足实时档预算")
    func highMuAdaptiveDomain() async throws {
        let module = PhotonStatisticsModule()
        // μ=50：σ=√(50×51)≈50.5 → n_max ≈ 50+15×50.5 ≈ 808
        #expect(PhotonStatsMath.nCompute(forMu: 50) >= 800)

        var values = ParamValues.defaults(for: module.params)
        values.sliders["mu"] = 50
        values.sliders["nshow"] = 60
        let result = try await module.compute(values, constants: try constants())
        guard case .bars(let bars) = result.charts[0] else {
            Issue.record("charts[0] 应为 bars"); return
        }
        #expect(bars.bars.count == 61, "nshow=60 → 61 柱")
        // ΣP 仍归一
        let sum = PhotonStatsMath.distribution(
            PhotonStatsMath.coherent, mu: 50, n: PhotonStatsMath.nCompute(forMu: 50)).reduce(0, +)
        #expect(abs(sum - 1) < 1e-9)

        try await expectComputeUnderBudget(
            module: module, values: values, constants: try constants(),
            budgetMillis: 16, runs: 11)
    }

    @Test("实时档预算：默认参数 compute 中位 < 16 ms")
    func computeBudget() async throws {
        let module = PhotonStatisticsModule()
        try await expectComputeUnderBudget(
            module: module,
            values: ParamValues.defaults(for: module.params),
            constants: try constants(),
            budgetMillis: 16)
    }
}
