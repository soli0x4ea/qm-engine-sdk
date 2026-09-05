import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W9-29 量子纠缠：纠缠度量与可分性。
/// Werner 解析曲线 + 负度/PPT 曲线 fixture 逐点对拍；Ginibre 系综按 RNG 政策做
/// 统计对拍（序列不复刻、分布等价——docs/RNG_SEED_POLICY.md）+ 物理律断言。
@Suite("W9 量子纠缠：纠缠度量与可分性")
struct EntanglementModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    // MARK: 计算核（物理律）

    @Test("物理律：Werner 阈值 p = 1/3（C = 0、min_eig = 0、N = 0）")
    func wernerThreshold() {
        let rho = EntanglementMath.wernerState(1.0 / 3.0)
        let pt = EntanglementMath.partialTransposeB(rho.re, rho.im)
        #expect(EntanglementMath.concurrence(re: rho.re, im: rho.im) < 1e-9)
        #expect(abs(EntanglementMath.minPptEig(ptRe: pt.re, ptIm: pt.im)) < 1e-6,
                "脚本断言口径：p = 1/3 处 PPT 最小本征值恰为 0")
        #expect(EntanglementMath.negativity(ptRe: pt.re, ptIm: pt.im) < 1e-9)
    }

    @Test("物理律：Wootters 数值并发度 = 解析 (3p−1)/2（p ∈ [1/3, 1] 六点）")
    func wernerAnalyticConcurrence() {
        for p in [1.0 / 3.0 + 1e-3, 0.4, 0.5, 0.6, 0.8, 1.0] {
            let rho = EntanglementMath.wernerState(p)
            let c = EntanglementMath.concurrence(re: rho.re, im: rho.im)
            let analytic = max(0.0, (3 * p - 1) / 2.0)
            #expect(abs(c - analytic) / max(analytic, 1e-12) < 1e-8,
                    "p = \(p): C = \(c) vs 解析 \(analytic)")
        }
        // p ≤ 1/3：C = 0
        for p in [0.0, 0.2, 1.0 / 3.0 - 1e-3] {
            let rho = EntanglementMath.wernerState(p)
            #expect(EntanglementMath.concurrence(re: rho.re, im: rho.im) < 1e-9)
        }
        // 单态极限 p = 1：C = 1、N = 1/2（‖ρ^T_B‖₁ = 2）
        let singlet = EntanglementMath.wernerState(1.0)
        let pt = EntanglementMath.partialTransposeB(singlet.re, singlet.im)
        #expect(abs(EntanglementMath.concurrence(re: singlet.re, im: singlet.im) - 1.0) < 1e-9)
        #expect(abs(EntanglementMath.negativity(ptRe: pt.re, ptIm: pt.im) - 0.5) < 1e-9)
    }

    @Test("物理律：纯态族 |ψ(θ)⟩ 并发度 = sin 2θ；纠缠熵峰值 1 ebit")
    func pureStateLaws() {
        func pureRho(_ theta: Double) -> (re: [Double], im: [Double]) {
            let c = cos(theta), s = sin(theta)
            let psi = [c, 0.0, 0.0, s]
            var re = [Double](repeating: 0, count: 16)
            for i in 0..<4 { for j in 0..<4 { re[i * 4 + j] = psi[i] * psi[j] } }
            return (re, EntanglementMath.zeros)
        }
        for theta in [Double.pi / 6, Double.pi / 4, Double.pi / 3] {
            let rho = pureRho(theta)
            let c = EntanglementMath.concurrence(re: rho.re, im: rho.im)
            // ρ^{1/2} 路径含幂迭代级数值噪声，容差 1e-8
            #expect(abs(c - sin(2 * theta)) < 1e-8, "θ = \(theta): C = \(c)")
        }
        #expect(abs(EntanglementMath.entEntropy(.pi / 4) - 1.0) < 1e-12)
        #expect(EntanglementMath.entEntropy(0) == 0)
        #expect(abs(EntanglementMath.entEntropy(.pi / 2)) < 1e-30)
    }

    @Test("物理律：Ginibre 样本 trace = 1 且半正定（前 300 样本不变量）")
    func ginibreInvariants() {
        for idx in 0..<300 {
            let rho = EntanglementMath.ginibreRho(
                EntanglementMath.makeNormalGenerator(seed: EntanglementModule.rngSeed,
                                                     sampleIndex: idx))
            let w = EntanglementMath.eigvals(rho.re, rho.im)
            let trace = w.reduce(0, +)
            #expect(abs(trace - 1.0) < 1e-12, "样本 \(idx): trace = \(trace)")
            #expect(w[0] > -1e-12, "样本 \(idx): min eig = \(w[0])（半正定）")
        }
    }

    /// 完整 4000 样本系综（测试内直算，与 compute 同种子同口径）。
    private func fullEnsemble() -> [(neg: Double, minPpt: Double, c: Double)] {
        (0..<EntanglementModule.fullSamples).map {
            EntanglementMath.ginibreSample(
                EntanglementMath.makeNormalGenerator(seed: EntanglementModule.rngSeed,
                                                     sampleIndex: $0))
        }
    }

    @Test("物理律：2×2 下 PPT ⇔ C = 0（4000 样本双向无一反例）")
    func pptImpliesSeparable() {
        let ensemble = fullEnsemble()
        var maxCPpt = 0.0
        var minCNpt = Double.infinity
        var nPpt = 0
        for s in ensemble {
            if s.minPpt >= -1e-10 {
                nPpt += 1
                maxCPpt = max(maxCPpt, s.c)
            } else {
                minCNpt = min(minCNpt, s.c)
            }
        }
        #expect(maxCPpt < 1e-6, "PPT 态最大 C = \(maxCPpt)（脚本断言口径）")
        #expect(minCNpt > 0, "NPT 态最小 C = \(minCNpt)")
        #expect(nPpt > 0 && nPpt < EntanglementModule.fullSamples)
    }

    // MARK: 统计对拍（RNG 政策：种子 = Python 原值，分布等价）

    @Test("统计对拍：PPT 占比 = Python 934/4000 ± 4.5σ")
    func pptFractionStatistics() {
        let ensemble = fullEnsemble()
        let nPpt = ensemble.filter { $0.minPpt >= -1e-10 }.count
        let pyCount = 934
        let q = Double(pyCount) / Double(EntanglementModule.fullSamples)
        let sigma = sqrt(Double(EntanglementModule.fullSamples) * q * (1 - q))
        #expect(abs(Double(nPpt - pyCount)) < 4.5 * sigma,
                "Swift PPT \(nPpt) vs Python \(pyCount)（σ = \(String(format: "%.1f", sigma))）")
    }

    @Test("统计对拍：NPT 组 ⟨N⟩/⟨C⟩ = fixture 散点均值（±0.05）")
    func scatterMeanStatistics() throws {
        let fixture = try ModuleFixture.load("量子纠缠_纠缠度量与可分性__v2022")
        let scatters = try fixture.scatters(2, 0)   // fig3：[PPT, NPT]
        #expect(scatters.count == 2)
        let nptFixture = scatters[1]
        let meanN = nptFixture.x.reduce(0, +) / Double(nptFixture.x.count)
        let meanC = nptFixture.y.reduce(0, +) / Double(nptFixture.y.count)

        let ensemble = fullEnsemble().filter { $0.minPpt < -1e-10 }
        let ourN = ensemble.map(\.neg).reduce(0, +) / Double(ensemble.count)
        let ourC = ensemble.map(\.c).reduce(0, +) / Double(ensemble.count)
        #expect(abs(ourN - meanN) < 0.05, "⟨N⟩: \(ourN) vs fixture \(meanN)")
        #expect(abs(ourC - meanC) < 0.05, "⟨C⟩: \(ourC) vs fixture \(meanC)")
    }

    // MARK: fixture 逐点对拍（Werner 201 点原生网格）

    @Test("fixture 对拍：Werner C/N/min_eig 三曲线 201 点（rel < 1e-8）")
    func wernerFixtureCurves() throws {
        let fixture = try ModuleFixture.load("量子纠缠_纠缠度量与可分性__v2022")
        let cFx = try fixture.line(0, 0, label: "Concurrence")
        let nFx = try fixture.line(0, 1, label: "Negativity")
        let minFx = try fixture.line(0, 1, label: "min eigenvalue")
        let ps = cFx.x   // 脚本口径 linspace(0,1,201)

        var ourC: [Double] = [], ourN: [Double] = [], ourMin: [Double] = []
        for p in ps {
            let rho = EntanglementMath.wernerState(p)
            let pt = EntanglementMath.partialTransposeB(rho.re, rho.im)
            ourC.append(EntanglementMath.concurrence(re: rho.re, im: rho.im))
            ourN.append(EntanglementMath.negativity(ptRe: pt.re, ptIm: pt.im))
            ourMin.append(EntanglementMath.minPptEig(ptRe: pt.re, ptIm: pt.im))
        }
        // C 曲线在 p ≤ 1/3 区域恒为 0（预期为 0 的点不能用纯相对判据）：
        // abs + rel 组合容差（fixture 存 9 位有效数字）
        func expectCombined(_ actual: [Double], _ expected: [Double], _ ctx: String) {
            #expect(actual.count == expected.count, "\(ctx): 点数不一致")
            for i in actual.indices {
                #expect(abs(actual[i] - expected[i]) <= 1e-8 + 1e-8 * abs(expected[i]),
                        "\(ctx) 第 \(i) 点：\(actual[i]) vs \(expected[i])")
            }
        }
        expectCombined(ourC, cFx.y, "Werner 并发度")
        expectCombined(ourN, nFx.y, "Werner 负度")
        expectCombined(ourMin, minFx.y, "Werner min eig")
    }

    // MARK: 模块契约

    @Test("模块契约：5 图表 + 帧栈 11 帧 4×4 + 散点两组 + 可复现（同种子逐位相同）")
    func moduleContractAndReproducibility() async throws {
        let module = EntanglementModule()
        var values = ParamValues.defaults(for: module.params)
        values.discretes["mode"] = "interactive"
        values.sliders["nInteractive"] = 200
        let c = try constants()

        let r1 = try await module.compute(values, constants: c)
        let r2 = try await module.compute(values, constants: c)

        #expect(r1.charts.count == 5)
        #expect(r1.charts[3] == r2.charts[3], "帧栈两次 compute 逐位相同（固定种子）")
        guard case .scatter(let s1) = r1.charts[4],
              case .scatter(let s2) = r2.charts[4] else {
            Issue.record("图 5 应为散点"); return
        }
        #expect(s1 == s2, "散点两次 compute 逐位相同（固定种子可复现）")
        #expect(s1.series.count == 2)

        guard case .frameStack(let fs) = r1.charts[3] else {
            Issue.record("图 4 应为帧栈"); return
        }
        #expect(fs.frames.count == 11)
        #expect(fs.frames.allSatisfy { $0.values.count == 4 && $0.values[0].count == 4 })
        #expect(fs.xTicks.count == 4 && fs.yTicks.count == 4)

        #expect(r1.summary.count == 5)
        #expect(r1.theory?.formulas.count == 4)
    }

    @Test("秒级档预算：交互 300 与完整 4000 均 < 2000 ms（强制口径）")
    func computeBudget() async throws {
        let module = EntanglementModule()
        let c = try constants()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: c, budgetMillis: 2000, runs: 3, warmup: 1)
        var full = ParamValues.defaults(for: module.params)
        full.discretes["mode"] = "full"
        try await expectComputeUnderBudget(
            module: module, values: full, constants: c,
            budgetMillis: 2000, runs: 3, warmup: 1)
    }
}
