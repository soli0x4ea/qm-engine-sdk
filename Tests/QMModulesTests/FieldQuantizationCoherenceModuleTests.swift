import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W6 笔记 35《光场量子化与相干态》：二阶相干度 g⁽²⁾(0) 三态对比 + 压缩真空 Wigner。
/// fixtures：光场量子化_二阶相干度__v2022.json（柱：相干 1.0 / 热态 1.99999986）。
@Suite("W6 二阶相干度")
struct FieldQuantizationCoherenceModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    // MARK: - 物理律（g⁽²⁾(0) 三态指纹）

    @Test("物理律：相干态 g⁽²⁾(0) = 1（泊松统计 / 散粒噪声极限），与光强 μ 无关")
    func coherentG2Law() {
        #expect(abs(SecondOrderCoherenceMath.g2Coherent(mu: 5) - 1) < 1e-9)
        #expect(abs(SecondOrderCoherenceMath.g2Coherent(mu: 0.5) - 1) < 1e-9)
        #expect(abs(SecondOrderCoherenceMath.g2Coherent(mu: 50) - 1) < 1e-9,
                "μ = 50：求和域自适应扩展（nmax = 808）后仍收敛到 1")
    }

    @Test("物理律：热态 g⁽²⁾(0) = 2（超泊松 / 聚束），与 μ 无关；三态排序 Fock < 相干 < 热")
    func thermalG2Law() {
        // 脚本口径 nmax=120：截断尾部 (5/6)¹²⁰ ≈ 3×10⁻¹⁰ 经矩比放大为 ~1×10⁻⁷ 偏差
        #expect(abs(SecondOrderCoherenceMath.g2Thermal(mu: 5) - 2) < 1e-6,
                "μ = 5（脚本 120 截断）：\(SecondOrderCoherenceMath.g2Thermal(mu: 5))")
        #expect(abs(SecondOrderCoherenceMath.g2Thermal(mu: 0.5) - 2) < 1e-9)
        // 扩展求和域（(5/6)⁶⁰⁰ ≈ 0）后收敛到 2——截断误差消失，极限律精确成立
        #expect(abs(SecondOrderCoherenceMath.g2Thermal(mu: 5, nmax: 600) - 2) < 1e-9,
                "扩展域收敛：\(SecondOrderCoherenceMath.g2Thermal(mu: 5, nmax: 600))")
        // μ = 50：截断尾部 (50/51)⁸⁰⁹ ≈ 1e-7 放大矩比误差至 ~3e-5
        #expect(abs(SecondOrderCoherenceMath.g2Thermal(mu: 50) - 2) < 1e-4)
        // 三态物理排序：反聚束 < 散粒噪声 < 聚束
        let fock = SecondOrderCoherenceMath.g2Fock(n: 1)
        let coh = SecondOrderCoherenceMath.g2Coherent(mu: 5)
        let th = SecondOrderCoherenceMath.g2Thermal(mu: 5)
        #expect(fock < coh && coh < th, "0 < 1 < 2：光场量子态的统计指纹")
    }

    @Test("物理律：Fock |n⟩ 的 g⁽²⁾(0) = 1 − 1/n——|1⟩ 反聚束极限为 0")
    func fockG2Law() {
        #expect(SecondOrderCoherenceMath.g2Fock(n: 1) == 0,
                "单光子态：Hanbury Brown-Twiss 零符合")
        #expect(abs(SecondOrderCoherenceMath.g2Fock(n: 2) - 0.5) < 1e-15)
        #expect(abs(SecondOrderCoherenceMath.g2Fock(n: 10) - 0.9) < 1e-15)
        for n in 1...50 {
            #expect(abs(SecondOrderCoherenceMath.g2Fock(n: n)
                        - (1.0 - 1.0 / Double(n))) < 1e-15)
            #expect(SecondOrderCoherenceMath.g2Fock(n: n) >= 0,
                    "概率比非负")
        }
    }

    // MARK: - 物理律（压缩真空 Wigner）

    @Test("物理律：Wigner 归一 ∫∫W dx dp = 1、恒正、峰值在原点 W(0,0) = 1/π")
    func wignerNormalizationLaw() {
        let r = 1.0
        #expect(abs(SecondOrderCoherenceMath.squeezedVacuumWigner(x: 0, p: 0, r: r)
                    - 1.0 / .pi) < 1e-15, "真空压缩不改变峰值 W(0,0) = 1/π")
        for (x, p) in [(-2.0, 0.3), (0.7, -1.4), (3.0, 3.0), (0, 0)] {
            #expect(SecondOrderCoherenceMath.squeezedVacuumWigner(x: x, p: p, r: r) >= 0,
                    "压缩真空是高斯态：Wigner 恒正")
            #expect(SecondOrderCoherenceMath.squeezedVacuumWigner(x: x, p: p, r: r)
                    <= 1.0 / .pi + 1e-15, "高斯峰值在原点")
        }
        // 模块同口径 320×320 网格上数值积分（±3σ 窗口 + 黎曼和，误差 < 1e-4）
        let n = 320
        let xMax = 3.0 * exp(-r), pMax = 3.0 * exp(r)
        let xs = Num.linspace(-xMax, xMax, count: n)
        let ps = Num.linspace(-pMax, pMax, count: n)
        var sum = 0.0
        for p in ps {
            for x in xs {
                sum += SecondOrderCoherenceMath.squeezedVacuumWigner(x: x, p: p, r: r)
            }
        }
        let integral = sum * (xs[1] - xs[0]) * (ps[1] - ps[0])
        #expect(abs(integral - 1) < 1e-4, "∫∫W = \(integral)")
    }

    @Test("物理律：1σ 椭圆 x²e^{2r} + p²e^{−2r} = 1 过 (e^{−r},0) 与 (0,e^{r})；最小不确定度")
    func squeezingEllipseLaw() {
        let r = 1.0
        let oneSigma = SecondOrderCoherenceMath.oneSigmaLevel()
        #expect(abs(oneSigma - exp(-1) / .pi) < 1e-15)
        #expect(abs(SecondOrderCoherenceMath.squeezedVacuumWigner(x: exp(-r), p: 0, r: r)
                    - oneSigma) < 1e-15, "长半轴端点 (e^{−r}, 0) 在 1σ 线上")
        #expect(abs(SecondOrderCoherenceMath.squeezedVacuumWigner(x: 0, p: exp(r), r: r)
                    - oneSigma) < 1e-15, "短半轴端点 (0, e^{r}) 在 1σ 线上")
        // 压缩 / 反压缩方差与海森伯积：⟨x²⟩⟨p²⟩ = 1/4（最小不确定度）
        let (varX, varP) = SecondOrderCoherenceMath.quadratureVariances(r: r)
        #expect(abs(varX - exp(-2 * r) / 2) < 1e-15)
        #expect(abs(varP - exp(2 * r) / 2) < 1e-15)
        #expect(abs(varX * varP - 0.25) < 1e-15, "压缩重新分配噪声而非违反不确定性")
        #expect(varX < 0.5 && varP > 0.5, "x 方向低于真空噪声 1/2（压缩）")
        // 压缩率守恒：⟨x²⟩·⟨p²⟩ 与 r 无关
        for r2 in [0.0, 0.5, 2.0] {
            let (vx, vp) = SecondOrderCoherenceMath.quadratureVariances(r: r2)
            #expect(abs(vx * vp - 0.25) < 1e-15)
        }
    }

    // MARK: - fixtures 对拍

    @Test("fixtures：μ = 5 柱——相干 1.0 / 热态 1.99999986（< 1e-8）")
    func barsMatchFixtures() throws {
        let fx = try ModuleFixture.load("光场量子化_二阶相干度__v2022")
        let bars = try fx.bars(0, 0)
        #expect(bars.count == 2, "脚本柱：相干态 + 热态（Fock 由理论卡文字给出）")
        #expect(bars[0].x == 0 && bars[1].x == 1)
        #expect(abs(bars[0].height - SecondOrderCoherenceMath.g2Coherent(mu: 5)) < 1e-8,
                "相干态 \(bars[0].height)")
        #expect(abs(bars[1].height - SecondOrderCoherenceMath.g2Thermal(mu: 5)) < 1e-8,
                "热态 \(bars[1].height)（120 截断的 1.99999986）")
        #expect(abs(bars[0].height - 1) < 1e-8 && abs(bars[1].height - 2) < 2e-7)
    }

    // MARK: - compute 结构与预算

    @Test("compute 结构：柱图 3 柱 + 320×320 等高线图 + 摘要 + 理论卡；秒级档预算 < 2 s")
    func computeStructureAndBudget() async throws {
        let module = SecondOrderCoherenceModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)

        guard case .bars(let bars) = result.charts[0] else {
            Issue.record("charts[0] 应为 bars"); return
        }
        #expect(bars.bars.count == 3)
        #expect(abs(bars.bars[0].value - 1) < 1e-8)
        #expect(abs(bars.bars[1].value - 2) < 2e-7, "热态 120 截断：\(bars.bars[1].value)")
        #expect(bars.bars[2].value == 0, "Fock |1⟩ 反聚束")

        guard case .contour(let contour) = result.charts[1] else {
            Issue.record("charts[1] 应为 contour"); return
        }
        #expect(contour.xGrid.count == 320 && contour.yGrid.count == 320)
        #expect(contour.values.count == 320 && contour.values[0].count == 320)
        #expect(contour.highlightLevels.count == 1, "1σ 高亮等值线")
        // 网格自检（偶数网格无原点格点；最近格点 |x| = a/319、|p| = b/319）：
        // 全局最大位于中心四格点，量级 1/π − 6×10⁻⁵（格点分辨率级偏差）
        let mid = 160
        let gridMax = contour.values.map { $0.max() ?? 0 }.max() ?? 0
        #expect(abs(gridMax - 1.0 / .pi) < 1e-3, "网格峰值 ≈ 1/π")
        #expect(contour.values[mid][mid] == gridMax, "峰值位于网格中心四格点")
        for row in contour.values {
            for w in row { #expect(w <= 1.0 / .pi + 1e-12 && w >= 0) }
        }

        #expect(result.summary.count == 5)
        #expect(result.theory?.formulas.count == 4)

        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 2000)
    }
}
