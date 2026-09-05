import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W9-46 矿物光学性质的第一性原理计算：Kramers–Kronig 反演。
/// Lorentz 振子闭式校验（Python 锚点）+ 六曲线 fixture 重采样对拍
/// + 求和律/吸收峰物理律。
@Suite("W9 矿物光学：Kramers-Kronig 反演")
struct KramersKronigModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    /// 脚本口径参数（fixture）：ω₀=5, ω_p=8, γ=0.6 eV。
    private var w0: Double { 5.0 }
    private var wp: Double { 8.0 }
    private var gam: Double { 0.6 }

    /// 脚本口径 KK 网格 + 反演（全量 12000²，供锚点测试复用）。
    private func runInversion() -> (w: [Double], eps2: [Double],
                                    eps1Exact: [Double], eps1Num: [Double]) {
        let w = Num.linspace(0.005, 40.0, count: KramersKronigModule.kkPoints)
        let eps2 = KramersKronigMath.lorentzEps2Grid(w, w0: w0, wp: wp, gamma: gam)
        let eps1Exact = KramersKronigMath.lorentzEps1Grid(w, w0: w0, wp: wp, gamma: gam)
        let eps1Num = KramersKronigMath.invertEpsilon1(w: w, eps2: eps2)
        return (w, eps2, eps1Exact, eps1Num)
    }

    // MARK: 物理律 / 脚本锚点

    @Test("脚本锚点：ε₁ max|Δ| ≈ 5.0993e-2；n(ω) max rel ≈ 1.3919e-3；ω=5 锚点")
    func scriptAnchors() {
        let r = runInversion()
        var maxAbsErr = 0.0
        var maxRelErrN = 0.0
        for i in r.w.indices {
            maxAbsErr = max(maxAbsErr, abs(r.eps1Num[i] - r.eps1Exact[i]))
            let ne = KramersKronigMath.refractiveIndex(eps1: r.eps1Exact[i], eps2: r.eps2[i]).n
            let nn = KramersKronigMath.refractiveIndex(eps1: r.eps1Num[i], eps2: r.eps2[i]).n
            maxRelErrN = max(maxRelErrN, abs(nn - ne) / (abs(ne) + 1e-9))
        }
        #expect(maxAbsErr < 6e-2 && maxAbsErr > 2e-2,
                "ε₁ max|Δ| = \(maxAbsErr)（Python 5.0993e-2，同算法应同量级）")
        #expect(maxRelErrN < 1.6e-3 && maxRelErrN > 1.0e-3,
                "n max rel = \(maxRelErrN)（Python 1.3919e-3）")

        // ω = 5 锚点：ε₁_exact = 0.8963、ε₁_num = 0.8992、ε₂ = 21.3266
        let idx = r.w.enumerated().min { abs($0.element - 5.0) < abs($1.element - 5.0) }!
            .offset
        #expect(abs(r.eps1Exact[idx] - 0.8963) < 5e-4, "ε₁_exact = \(r.eps1Exact[idx])")
        #expect(abs(r.eps1Num[idx] - 0.8992) < 5e-4, "ε₁_num = \(r.eps1Num[idx])")
        #expect(abs(r.eps2[idx] - 21.3266) < 5e-4, "ε₂ = \(r.eps2[idx])")
    }

    @Test("物理律：求和律 ε₁(0⁺) = 1 + ω_p²/ω₀² = 3.56；ε₂ 峰位落在 [ω₀−γ, ω₀]")
    func sumRuleAndPeak() {
        // 低频极限：解析 ε₁(0.005) ≈ 1 + ω_p²/ω₀²（ω → 0 时 KK 的零频行为）
        let eps1Low = KramersKronigMath.lorentzEps1(0.005, w0: w0, wp: wp, gamma: gam)
        #expect(abs(eps1Low - (1.0 + wp * wp / (w0 * w0))) < 1e-3,
                "ε₁(0⁺) = \(eps1Low) vs 求和律 3.56")

        // 吸收峰：ε₂ 最大值位置在共振 ω₀ 附近（阻尼使峰位略低于 ω₀）
        let wFine = Num.linspace(3.0, 7.0, count: 4000)
        let eps2 = wFine.map { KramersKronigMath.lorentzEps2($0, w0: w0, wp: wp, gamma: gam) }
        let peakW = wFine[eps2.enumerated().max { $0.element < $1.element }!.offset]
        #expect(peakW > w0 - gam && peakW <= w0 + 1e-9,
                "ε₂ 峰位 \(peakW) 应落在 [ω₀−γ, ω₀]")

        // n² − κ² = ε₁、2nκ = ε₂（复折射率定义自洽）
        let ri = KramersKronigMath.refractiveIndex(
            eps1: KramersKronigMath.lorentzEps1(5.0, w0: w0, wp: wp, gamma: gam),
            eps2: KramersKronigMath.lorentzEps2(5.0, w0: w0, wp: wp, gamma: gam))
        let eps1 = KramersKronigMath.lorentzEps1(5.0, w0: w0, wp: wp, gamma: gam)
        let eps2At5 = KramersKronigMath.lorentzEps2(5.0, w0: w0, wp: wp, gamma: gam)
        #expect(abs(ri.n * ri.n - ri.kappa * ri.kappa - eps1) < 1e-12)
        #expect(abs(2.0 * ri.n * ri.kappa - eps2At5) < 1e-12)
    }

    // MARK: fixture 重采样对拍（512 点管线产物 → 模块绘图网格）

    @Test("fixture 对拍：ε₂/ε₁解析/ε₁KK + n/n/κ 六曲线 512 点（rel < 5e-3）")
    func fixtureCurves() throws {
        let fixture = try ModuleFixture.load("46_矿物光学性质的第一性原理计算_KramersKronig__v2022")
        let eps2Fx = try fixture.line(0, 0, label: "input")
        let eps1AFx = try fixture.line(0, 0, label: "analytic")
        let eps1NFx = try fixture.line(0, 0, label: "from KK inversion")
        let nAFx = try fixture.line(0, 1, label: "analytic")
        let nNFx = try fixture.line(0, 1, label: "from KK + ")
        let kFx = try fixture.line(0, 1, label: "kappa")

        let w = eps2Fx.x   // fixture 绘图网格（0.05…12 eV）
        let eps2 = w.map { KramersKronigMath.lorentzEps2($0, w0: w0, wp: wp, gamma: gam) }
        let eps1A = w.map { KramersKronigMath.lorentzEps1($0, w0: w0, wp: wp, gamma: gam) }
        let inv = runInversion()
        let eps1N = w.map { KramersKronigMath.interp($0, xs: inv.w, ys: inv.eps1Num) }
        let nA = zip(eps1A, eps2).map { KramersKronigMath.refractiveIndex(eps1: $0, eps2: $1).n }
        let nN = zip(eps1N, eps2).map { KramersKronigMath.refractiveIndex(eps1: $0, eps2: $1).n }
        let kap = zip(eps1A, eps2).map { KramersKronigMath.refractiveIndex(eps1: $0, eps2: $1).kappa }

        // 512 点为管线重采样，取相对容差 5e-3（重采样 + 插值误差主导）
        expectResampled(moduleX: w, moduleY: eps2, fx: eps2Fx.x, fy: eps2Fx.y,
                        tolerance: 5e-3, "ε₂ 输入")
        expectResampled(moduleX: w, moduleY: eps1A, fx: eps1AFx.x, fy: eps1AFx.y,
                        tolerance: 5e-3, "ε₁ 解析")
        expectResampled(moduleX: w, moduleY: eps1N, fx: eps1NFx.x, fy: eps1NFx.y,
                        tolerance: 5e-3, "ε₁ KK 反演")
        expectResampled(moduleX: w, moduleY: nA, fx: nAFx.x, fy: nAFx.y,
                        tolerance: 5e-3, "n 解析")
        expectResampled(moduleX: w, moduleY: nN, fx: nNFx.x, fy: nNFx.y,
                        tolerance: 5e-3, "n KK + ε₂")
        expectResampled(moduleX: w, moduleY: kap, fx: kFx.x, fy: kFx.y,
                        tolerance: 5e-3, "κ(ω)")
    }

    // MARK: 模块契约

    @Test("模块契约：双图 3+3 序列各 600 点 + 摘要 5 卡")
    func moduleContract() async throws {
        let module = KramersKronigModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())

        #expect(result.charts.count == 2)
        guard case .lineSeries(let c1) = result.charts[0],
              case .lineSeries(let c2) = result.charts[1] else {
            Issue.record("两图均应为曲线"); return
        }
        #expect(c1.series.count == 3 && c2.series.count == 3)
        #expect(c1.series.allSatisfy { $0.points.count == 600 })
        #expect(c2.series.allSatisfy { $0.points.count == 600 })

        #expect(result.summary.count == 5)
        #expect(result.summary[4].value.contains("1 + ω_p²/ω₀²"))
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("秒级档预算：compute < 2000 ms（12000² vDSP 反演，实测登记）")
    func computeBudget() async throws {
        let module = KramersKronigModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 2000, runs: 3, warmup: 1)
    }
}
