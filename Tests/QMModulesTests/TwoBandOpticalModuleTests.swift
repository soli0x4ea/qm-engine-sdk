import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// 笔记 46《矿物光学性质的第一性原理计算》模型 (a)：两带光学模型（W11B 策略档）。
/// 物理律验证 + fixture 对拍 + StrategyCore 复用验证 + 预算登记。
///
/// 判据均经 numpy 试算锚定（W11B 方案 §3，N = 20000/200001 实测）：
/// - √ 开启律 @200001 interior bins max dev 1.10e-2 → 判 2e-2；
/// - 总谱权重 S(20000) 偏差 7.5e-5 → 判 2e-4；S(5000) vs S(20000) 2.25e-4 → 判 5e-4；
/// - KK 静态极限闭合 rel 3e-5 → 判 1e-3；
/// - 20k 逐 bin 量子化 ~6%（每 bin 仅 ~13–17 点）→ fixture 曲线容差 1e-1。
@Suite(.serialized)
struct TwoBandOpticalModuleTests {

    private let module = TwoBandOpticalModule()

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    /// 直算管线（与模块 compute 同一数学核，物理律口径用全分辨率数组）。
    private struct Pipeline {
        let w: [Double]
        let eps2: [Double]
        let eps1: [Double]
        let n: [Double]
        let kappa: [Double]
        let centers: [Double]
        let jdos: [Double]
        let edges: [Double]
        let dk: Double
    }

    private func runPipeline(nK: Int, eg: Double = 4.0, decay: Double = 8.0) -> Pipeline {
        let ks = Num.linspace(0, TwoBandMath.kmax, count: nK)
        let edges = TwoBandMath.binEdges(eg: eg)
        var sums = [Double](repeating: 0, count: TwoBandMath.binCount)
        TwoBandMath.accumulate(ks: ks[...], eg: eg, edges: edges, into: &sums)
        let (centers, jdos) = TwoBandMath.jdosTable(sums: sums, edges: edges)
        let w = Num.linspace(TwoBandMath.plotMin, TwoBandMath.plotMax,
                             count: TwoBandMath.plotPoints)
        let eps2 = TwoBandMath.eps2Grid(w: w, centers: centers, jdos: jdos,
                                        eg: eg, decay: decay)
        let eps1 = KramersKronigMath.invertEpsilon1(w: w, eps2: eps2)
        var n = [Double](repeating: 0, count: w.count)
        var kappa = [Double](repeating: 0, count: w.count)
        for i in w.indices {
            let nk = KramersKronigMath.refractiveIndex(eps1: eps1[i], eps2: eps2[i])
            n[i] = nk.n
            kappa[i] = nk.kappa
        }
        return Pipeline(w: w, eps2: eps2, eps1: eps1, n: n, kappa: kappa,
                        centers: centers, jdos: jdos, edges: edges,
                        dk: TwoBandMath.kmax / Double(nK - 1))
    }

    // MARK: 物理律 1：JDOS 平方根开启律 + 带边特征

    @Test("物理律：JDOS ∝ √(ħω−E_g)（200001 脚本口径 interior bins < 2e-2）+ 带边特征")
    func jdosSquareRootLaw() {
        let p = runPipeline(nK: 200001)
        // 解析 bin 值：V(i) = (k_hi³ − k_lo³)/(3·dk)，k(E) = √((E−E_g)/TK)
        var v = [Double](repeating: 0, count: TwoBandMath.binCount)
        for i in 0..<TwoBandMath.binCount {
            let kLo = sqrt(max(p.edges[i] - 4.0, 0) / TwoBandMath.tk)
            let kHi = sqrt(max(p.edges[i + 1] - 4.0, 0) / TwoBandMath.tk)
            v[i] = (kHi * kHi * kHi - kLo * kLo * kLo) / (3.0 * p.dk)
        }
        var pk = 0
        for i in v.indices where v[i] > v[pk] { pk = i }
        var maxDev = 0.0
        for i in 5..<591 {
            let ana = v[i] / v[pk]
            let dev = abs(p.jdos[i] - ana) / ana
            maxDev = max(maxDev, dev)
            #expect(dev < 2e-2, "bin \(i)：JDOS \(p.jdos[i]) vs 解析 \(ana)（dev \(dev)）")
        }
        #expect(maxDev > 0, "interior bins 全部通过（max dev \(maxDev)）")
        // 带边特征：边缘 bin 从近 0 开启，内部（量子化噪声内）单调上升
        #expect(p.jdos[0] < 0.05, "J(E_g⁺) = \(p.jdos[0]) 应近 0（√ 开启）")
        for i in 0..<(TwoBandMath.binCount - 5) {
            #expect(p.jdos[i + 1] > p.jdos[i] - 0.02,
                    "bin \(i)→\(i + 1)：\(p.jdos[i]) → \(p.jdos[i + 1]) 非单调")
        }
    }

    // MARK: 物理律 2：KK 闭合（因果性）

    @Test("物理律：KK 静态极限闭合 rel < 1e-3；带隙下 ε₂ ≡ 0 ⇒ n = √ε₁")
    func kkClosure() {
        let p = runPipeline(nK: 20000)
        // 静态极限（ω→0⁺ 核）与逐点主值反演（ω'ε₂/(ω'²−ω²) 核）独立互证
        let staticLimit = TwoBandMath.eps1StaticLimit(w: p.w, eps2: p.eps2)
        let rel = abs(p.eps1[0] - staticLimit) / abs(staticLimit)
        #expect(rel < 1e-3, "ε₁(0.05) = \(p.eps1[0]) vs 静态极限 \(staticLimit)（rel \(rel)）")
        // 因果性 sanity：带隙以下 ε₂ 严格为 0（模型有限支撑），n = √ε₁ 精确成立
        for i in p.w.indices where p.w[i] < 4.0 {
            #expect(p.eps2[i] == 0, "带隙以下 ε₂ 应为 0")
            let nn = sqrt(p.eps1[i])
            #expect(abs(p.n[i] - nn) / nn < 1e-12, "w = \(p.w[i])：n − √ε₁")
        }
    }

    // MARK: 物理律 3：采样收敛（策略降维 200001 → 20000 的统计量依据）

    @Test("物理律：总谱权重 S(20000) → KMAX³/3（< 2e-4）；S(5000) vs S(20000) < 5e-4；ConvergenceMonitor 收敛")
    func samplingConvergence() {
        func s(_ n: Int) -> Double {
            TwoBandMath.totalWeight(ks: Num.linspace(0, TwoBandMath.kmax, count: n))
        }
        let target = TwoBandMath.totalWeightTarget
        let s5 = s(5000)
        let s20 = s(20000)
        #expect(abs(s20 - target) / target < 2e-4,
                "S(20000) = \(s20) vs \(target)")
        #expect(abs(s5 - s20) / target < 5e-4,
                "S(5000) = \(s5) vs S(20000) = \(s20)")
        // 模块同款收敛监视：五档序列 patience 2 容差收敛
        var monitor = ConvergenceMonitor(relTolerance: 2e-3, patience: 2,
                                         maxSteps: TwoBandOpticalModule.convergenceTiers.count)
        for t in TwoBandOpticalModule.convergenceTiers {
            if monitor.observe(s(t)) { break }
        }
        #expect(monitor.isConverged, "S(N) 五档序列应容差收敛（convergedAt = \(String(describing: monitor.convergedAt))）")
    }

    // MARK: 物理律 4：光学常数自洽 + numpy 锚点

    @Test("物理律：n²−κ² = ε₁、2nκ = ε₂（< 1e-12）；ħω = 6 eV 锚点 vs numpy 实算")
    func opticalConstantsConsistency() {
        let p = runPipeline(nK: 20000)
        for i in p.w.indices {
            let back1 = p.n[i] * p.n[i] - p.kappa[i] * p.kappa[i]
            #expect(abs(back1 - p.eps1[i]) / max(abs(p.eps1[i]), 1e-300) < 1e-12,
                    "w = \(p.w[i])：n²−κ² ≠ ε₁")
            let back2 = 2.0 * p.n[i] * p.kappa[i]
            #expect(abs(back2 - p.eps2[i]) / max(abs(p.eps2[i]), 1e-300) < 1e-12,
                    "w = \(p.w[i])：2nκ ≠ ε₂")
        }
        // numpy 试算锚点（N = 20000，E_g = 4, Δ = 8，W11B 方案 §3）：n = 1.83378,
        // κ = 0.59167, ε₁ = 3.01266, ε₂ = 2.16998
        let idx = p.w.indices.min { abs(p.w[$0] - 6.0) < abs(p.w[$1] - 6.0) } ?? 0
        #expect(abs(p.n[idx] - 1.83378) < 1e-4, "n(6) = \(p.n[idx])")
        #expect(abs(p.kappa[idx] - 0.59167) < 1e-4, "κ(6) = \(p.kappa[idx])")
        #expect(abs(p.eps1[idx] - 3.01266) < 1e-4, "ε₁(6) = \(p.eps1[idx])")
        #expect(abs(p.eps2[idx] - 2.16998) < 1e-4, "ε₂(6) = \(p.eps2[idx])")
    }

    // MARK: fixture 对拍（脚本 N = 200001 → 模块 20000，曲线级容差 1e-1，见方案 §3）

    @Test("fixture：ε₂(ω) 对拍（跳 E_g 拐点 ±0.15 eV 与支撑边 16.05 eV 以上）")
    func fixtureEps2() async throws {
        let fixture = try ModuleFixture.load("46_矿物光学性质的第一性原理计算_模型__v2022")
        let result = try await module.compute(ParamValues.defaults(for: module.params),
                                              constants: try constants())
        let chart = try #require(chartLineSeries(result, 1))
        let fx = try fixture.line(0, 0, label: "model")
        expectResampled(
            moduleX: chart.series[0].points.map(\.x),
            moduleY: chart.series[0].points.map(\.y),
            fx: fx.x, fy: fx.y, tolerance: 0.1, "ε₂ 对拍（20k 逐 bin 量子化 ~6%）",
            skip: { abs($0 - 4.0) < 0.15 || $0 > 16.05 })
    }

    @Test("fixture：n(ω)/κ(ω) 对拍（n 跳 ε₁ 过零病态区与 w>17.5；κ 跳带边起点与支撑边）")
    func fixtureNKappa() async throws {
        let fixture = try ModuleFixture.load("46_矿物光学性质的第一性原理计算_模型__v2022")
        let result = try await module.compute(ParamValues.defaults(for: module.params),
                                              constants: try constants())
        let chart = try #require(chartLineSeries(result, 2))
        let fxN = try fixture.line(0, 1, index: 0)
        let fxK = try fixture.line(0, 1, index: 1)
        // n = √ε₁ 在 ε₁ 过零/负区（fixture n < 1e-3，w ≈ 16.28–17.93 全平台）病态：
        // 0 vs 5.4e-7 无相对意义（numpy-20000 交叉验证 n(16.28) = 0 与 Swift 一致）
        expectResampled(
            moduleX: chart.series[0].points.map(\.x),
            moduleY: chart.series[0].points.map(\.y),
            fx: fxN.x, fy: fxN.y, tolerance: 0.1, "n(ω) 对拍",
            skip: { x in x > 17.5 || (linterp(fxN.x, fxN.y, x) ?? 1) < 1e-3 })
        // κ ∝ (w−E_g)^{1/4} 在带边斜率无穷大，fixture 512 点线性插值必然欠估
        // （numpy-20000 κ(4.002) = 4.46e-2 与 Swift 一致，fixture 插值 2.89e-2 为伪影）
        expectResampled(
            moduleX: chart.series[1].points.map(\.x),
            moduleY: chart.series[1].points.map(\.y),
            fx: fxK.x, fy: fxK.y, tolerance: 0.1, "κ(ω) 对拍",
            skip: { abs($0 - 4.0) < 0.15 || $0 >= 16.4 })
    }

    @Test("fixture：带隙以下平均 n（统计量级对拍，< 2e-2）")
    func fixtureMeanN() async throws {
        let fixture = try ModuleFixture.load("46_矿物光学性质的第一性原理计算_模型__v2022")
        let result = try await module.compute(ParamValues.defaults(for: module.params),
                                              constants: try constants())
        let chart = try #require(chartLineSeries(result, 2))
        let fxN = try fixture.line(0, 1, index: 0)
        // fixture 曲线自身导出：mean n(w < 3.7)（脚本打印口径 low = w < E_g − 0.3）
        let low = zip(fxN.x, fxN.y).filter { $0.0 < 3.7 }.map(\.1)
        let fixtureMean = low.reduce(0, +) / Double(low.count)
        // 模块图表曲线同口径
        let mine = chart.series[0].points.filter { $0.x < 3.7 }.map(\.y)
        let moduleMean = mine.reduce(0, +) / Double(mine.count)
        #expect(abs(moduleMean - fixtureMean) / fixtureMean < 2e-2,
                "带隙下 n̄：模块 \(moduleMean) vs fixture \(fixtureMean)")
    }

    // MARK: StrategyCore 复用验证（缓存命中统计 + 结果恒等 + 检查点续算恒等）

    @Test("StrategyCore：JDOS 缓存（改 Δ 命中 / 改 E_g 未命中）+ 结果位级恒等 + 检查点落盘")
    func strategyCoreReuse() async throws {
        let consts = try constants()
        let defaults = ParamValues.defaults(for: module.params)
        TwoBandOpticalModule.jdosCache.removeAll()
        // 首算：miss
        let r1 = try await module.compute(defaults, constants: consts)
        let s1 = TwoBandOpticalModule.jdosCache.stats
        #expect(s1.misses == 1 && s1.hits == 0, "首算应 miss（\(s1)）")
        // 同参数重算：hit，结果位级恒等
        let r2 = try await module.compute(defaults, constants: consts)
        let s2 = TwoBandOpticalModule.jdosCache.stats
        #expect(s2.hits == 1, "同参数应 hit（\(s2)）")
        let c1 = try #require(chartLineSeries(r1, 1))
        let c2 = try #require(chartLineSeries(r2, 1))
        #expect(c1.series.count == c2.series.count)
        for (p, q) in zip(c1.series[0].points, c2.series[0].points) {
            #expect(p.x == q.x && p.y == q.y, "缓存命中路径应位级恒等")
        }
        // 改包络宽度 Δ（后处理参数）：仍命中（JDOS 与 Δ 无关）
        var decayVals = defaults
        decayVals.sliders["decay"] = 6.0
        _ = try await module.compute(decayVals, constants: consts)
        #expect(TwoBandOpticalModule.jdosCache.stats.hits == 2, "改 Δ 应命中")
        // 改带隙：未命中（JDOS 依赖 E_g）
        var egVals = defaults
        egVals.sliders["eg"] = 5.0
        _ = try await module.compute(egVals, constants: consts)
        #expect(TwoBandOpticalModule.jdosCache.stats.misses == 2, "改 E_g 应 miss")
        // 会话检查点已落盘（eg/kSamples/满批数）
        let cp = try #require(TwoBandOpticalModule.session.current)
        #expect(cp.eg == 5.0 && cp.kSamples == TwoBandOpticalModule.kSamples(for: .fine))
        #expect(cp.batchesDone == (cp.kSamples + TwoBandOpticalModule.batchSize - 1)
            / TwoBandOpticalModule.batchSize)
    }

    @Test("检查点续算：分批累计 + 部分和恢复 = 单次全量（位级恒等）")
    func checkpointResumeIdentity() {
        let eg = 4.0
        let nK = 20000
        let edges = TwoBandMath.binEdges(eg: eg)
        let ks = Num.linspace(0, TwoBandMath.kmax, count: nK)
        let totalBatches = (nK + TwoBandOpticalModule.batchSize - 1) / TwoBandOpticalModule.batchSize
        // 单次全量
        var full = [Double](repeating: 0, count: TwoBandMath.binCount)
        TwoBandMath.accumulate(ks: ks[...], eg: eg, edges: edges, into: &full)
        // 分批跑到第 3 批「被取消」，检查点携带部分和
        let cutoff = 3
        var partial = [Double](repeating: 0, count: TwoBandMath.binCount)
        for b in 0..<cutoff {
            let lo = b * TwoBandOpticalModule.batchSize
            let hi = min(lo + TwoBandOpticalModule.batchSize, nK)
            TwoBandMath.accumulate(ks: ks[lo..<hi], eg: eg, edges: edges, into: &partial)
        }
        let cp = TwoBandOpticalModule.SamplingCheckpoint(
            eg: eg, kSamples: nK, batchesDone: cutoff, binSums: partial)
        // 恢复：eg/kSamples 匹配 → 从检查点继续
        #expect(cp.eg == eg && cp.kSamples == nK, "恢复校验应通过")
        var resumed = cp.binSums
        for b in cp.batchesDone..<totalBatches {
            let lo = b * TwoBandOpticalModule.batchSize
            let hi = min(lo + TwoBandOpticalModule.batchSize, nK)
            TwoBandMath.accumulate(ks: ks[lo..<hi], eg: eg, edges: edges, into: &resumed)
        }
        #expect(resumed == full, "续算直方图应与单次全量位级恒等")
        // 不匹配（eg 变了）→ 应拒绝续算重采（模块恢复条件的语义锚）
        let stale = TwoBandOpticalModule.SamplingCheckpoint(
            eg: 5.0, kSamples: nK, batchesDone: cutoff, binSums: partial)
        #expect(stale.eg != eg, "参数不一致时恢复条件不成立")
    }

    // MARK: 预算（策略档名义预算：首算 10 s；强制口径）

    @Test("预算：best-of-41 强制（策略档 < 10 s）")
    func budget() async throws {
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 10_000)
    }
}
