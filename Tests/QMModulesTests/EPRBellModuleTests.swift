import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W5-11 EPR-Bell：CHSH / LHV 界 / Tsirelson 界（脚本六断言复现 +
/// fixtures 图 0 双轴 512 点曲线对拍 + 20 万定种子采样统计）。秒级档。
@Suite("W5 EPR 与 Bell 不等式")
struct EPRBellModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    @Test("LHV 枚举：16 种 ±1 取值组合 max|S| = 2（脚本 [1]，1e-12）")
    func lhvBound() {
        #expect(abs(EPRBellMath.lhvMaxS() - 2.0) < 1e-12)
    }

    @Test("自旋单态：S(π/4) = −2√2 解析值；1000 点网格极值差 < 1e-3（脚本 [2]）")
    func singlet() {
        #expect(abs(EPRBellMath.singletS(.pi / 4) + 2 * sqrt(2)) < 1e-12)
        let grid = Num.linspace(0, .pi / 2, count: 1000)
        var best = 0.0
        for phi in grid { best = max(best, abs(EPRBellMath.singletS(phi))) }
        #expect(abs(best - 2 * sqrt(2)) < 1e-3, "|S|max = \(best)")
        // S(0) = −2 恰落 LHV 界上
        #expect(abs(EPRBellMath.singletS(0) + 2) < 1e-12)
    }

    @Test("光子偏振：S(π/8) = +2√2（Bell 角 22.5°）；2000 点网格极值差 < 1e-3（脚本 [3]）")
    func photon() {
        #expect(abs(EPRBellMath.photonS(.pi / 8) - 2 * sqrt(2)) < 1e-12)
        let grid = Num.linspace(0, .pi / 4, count: 2000)
        var best = 0.0
        for t in grid { best = max(best, EPRBellMath.photonS(t)) }
        #expect(abs(best - 2 * sqrt(2)) < 1e-3, "S(θ*) = \(best)")
        #expect(abs(EPRBellMath.photonS(0) - 2) < 1e-12)
    }

    @Test("Bell 角配置 (0°,45°|22.5°,67.5°) → S = 2√2")
    func bellAngles() {
        let d2r = Double.pi / 180
        let s = EPRBellMath.chsh(
            EPRBellMath.E_phipol(0, 22.5 * d2r),
            EPRBellMath.E_phipol(0, 67.5 * d2r),
            EPRBellMath.E_phipol(45 * d2r, 22.5 * d2r),
            EPRBellMath.E_phipol(45 * d2r, 67.5 * d2r))
        #expect(abs(s - 2 * sqrt(2)) < 1e-12, "S = \(s)")
    }

    @Test("Tsirelson 采样（20 万、种子 42）：无一越界，极值充分逼近界（脚本 [4]）")
    func tsirelsonSampling() {
        let n = 200_000
        let s = EPRBellMath.tsirelsonSampling(n: n)
        #expect(s.violations == 0, "越界样本数 = \(s.violations)")
        #expect(s.maxS <= EPRBellMath.tsirelson + 1e-9)
        #expect(s.minS >= -EPRBellMath.tsirelson - 1e-9)
        #expect(s.maxS > 2.7, "maxS = \(s.maxS)，20 万采样应充分逼近 2√2")
        #expect(s.minS < -2.7, "minS = \(s.minS)")
        // 定种子可复现性：同参数两次调用逐位一致（RNG 规范 docs/RNG_SEED_POLICY.md）
        let s2 = EPRBellMath.tsirelsonSampling(n: n)
        #expect(s2.maxS == s.maxS && s2.minS == s.minS && s2.violations == s.violations)
    }

    @Test("Bell (1964)：lhs = √2/2 > rhs = 1 − √2/2（脚本 [4b]）")
    func bell1964() {
        let b = EPRBellMath.bell1964()
        #expect(b.lhs > b.rhs)
        #expect(abs(b.lhs - sqrt(2) / 2) < 1e-12)
        #expect(abs(b.rhs - (1 - sqrt(2) / 2)) < 1e-12)
    }

    @Test("Werner 阈值：p* ≈ 1/√2（网格步长 1e-4 内，脚本 [4c]）")
    func werner() {
        let w = EPRBellMath.wernerThreshold()
        #expect(w.p.count == 10001)
        #expect(abs(w.pThr - 1 / sqrt(2)) < 1e-3, "p* = \(w.pThr)")
        // 阈值前的最后一个点不违反、阈值点开始违反
        let idx = w.p.firstIndex(of: w.pThr) ?? 0
        #expect(w.s[idx] > 2 && w.s[idx - 1] <= 2)
    }

    @Test("fixture 曲线对拍：自旋单态 + 光子偏振（图 0 双轴，各 512 点）")
    func curveMatchFixture() throws {
        let fx = try ModuleFixture.load("EPR与Bell不等式_CHSH__v2022")
        let d2r = Double.pi / 180
        // x 为 9 位有效数字舍入的度值：解析公式复算，绝对容差 1e-8（S 量级 ~1）
        let (x1, y1) = try fx.line(0, 0, index: 0)
        #expect(x1.count == 512)
        #expect(abs(x1.first! - 0) < 1e-12 && abs(x1.last! - 90) < 1e-9)
        for i in x1.indices {
            let s = EPRBellMath.singletS(x1[i] * d2r)
            #expect(abs(s - y1[i]) < 1e-8, "singlet 第 \(i) 点：\(s) vs \(y1[i])")
        }
        let (x2, y2) = try fx.line(0, 1, index: 0)
        #expect(x2.count == 512)
        #expect(abs(x2.first! - 0) < 1e-12 && abs(x2.last! - 45) < 1e-9)
        // 光子关联含 cos² 高阶导数，x 9 位舍入放大 → 绝对容差 1e-7（S 量级 ~1）
        for i in x2.indices {
            let s = EPRBellMath.photonS(x2[i] * d2r)
            #expect(abs(s - y2[i]) < 1e-7, "photon 第 \(i) 点：\(s) vs \(y2[i])")
        }
    }

    @Test("compute 输出结构：3 图（2 lineSeries 500 点 + 1 scatter）+ 摘要 7 条 + 理论卡")
    func computeStructure() async throws {
        let module = EPRBellModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 3)
        guard case .lineSeries(let c0) = result.charts[0],
              case .lineSeries(let c1) = result.charts[1],
              case .scatter(let c2) = result.charts[2] else {
            Issue.record("应为 2 lineSeries + 1 scatter"); return
        }
        // 单态 1000 点 stride 2 → 500；光子 2000 点 stride 4 → 500
        #expect(c0.series[0].points.count == 500)
        #expect(c1.series[0].points.count == 500)
        #expect(c0.referenceLines.count == 5, "LHV±2 + Tsirelson±2√2 + φ*")
        #expect(c1.referenceLines.count == 3, "LHV 2 + Tsirelson 2√2 + θ*")
        // 实验散点：3 个实验、LHV/Tsirelson 竖参考线
        #expect(c2.series.count == 3)
        #expect(c2.series.map(\.name) == EPRBellExperiments.all.map(\.label))
        let expS = c2.series.map { $0.points.first!.x }
        #expect(expS == [2.697, 2.25, 2.42])
        #expect(c2.referenceLines.count == 2)
        #expect(result.summary.count == 7)
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("秒级档预算：compute 中位 < 2000 ms（默认 N = 20 万）")
    func computeBudget() async throws {
        let module = EPRBellModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 2000)
    }
}
