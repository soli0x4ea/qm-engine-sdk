import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W5-12 晶体电子 Kronig-Penney：能带方程扫描 + 一维态密度（Python 端
/// 20000 点扫描重算对拍 fixtures 图 1 g(E) 512 点；无 [check] 断言，公式复算口径）。
@Suite("W5 晶体电子：Kronig-Penney 能带")
struct KronigPenneyModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    /// 与模块/脚本同口径的 20000 点扫描（P=1, a=1, Emax=160）。
    private func scan(P: Double = 1, emax: Double = 160)
        -> (E: [Double], allowed: [Bool], kRed: [Double]) {
        let E = Num.linspace(1e-4, emax, count: 20000)
        var allowed = [Bool](repeating: false, count: E.count)
        var kRed = [Double](repeating: 0, count: E.count)
        for i in E.indices {
            let f = KronigPenneyMath.f(E[i], P: P, a: 1)
            allowed[i] = abs(f) <= 1
            kRed[i] = acos(min(max(f, -1), 1))
        }
        return (E, allowed, kRed)
    }

    @Test("f(E) 解析核：xa = π/2 处 f = 2/π（P=1）；E→0⁺ 极限 f → 1 + P")
    func fKernel() {
        let E = KronigPenneyMath.hbar2Over2m * pow(.pi / 2, 2)  // xa = π/2, a = 1
        #expect(abs(KronigPenneyMath.f(E, P: 1, a: 1) - 2 / .pi) < 1e-12)
        // xa → 0：cos < 1、sin(xa)/(xa) < 1（单调下方趋近）⇒ f ↑ 1+P = 2（禁带）
        let f0 = KronigPenneyMath.f(1e-4, P: 1, a: 1)
        #expect(f0 < 2 && abs(f0 - 2) < 1e-3, "f(0⁺) → 2，实测 \(f0)")
        // 深带底（高 E）：f 振荡项主导 → 允许带与禁带交替
        #expect(KronigPenneyMath.kReduced(E, P: 1, a: 1) >= 0)
        #expect(KronigPenneyMath.kReduced(E, P: 1, a: 1) <= .pi + 1e-12)
    }

    @Test("P=1, Emax=160：允许带恰 2 段 [6.504, 37.602] / [51.411, 150.408]，第一禁带 13.809 eV")
    func bandSegments() {
        let (E, allowed, _) = scan()
        let segs = KronigPenneyMath.segments(allowed: allowed, E: E)
        #expect(segs.count == 2, "段数 = \(segs.count)")
        #expect(abs(segs[0].lo - 6.5044) < 5e-3)
        #expect(abs(segs[0].hi - 37.602) < 5e-3)
        #expect(abs(segs[1].lo - 51.4106) < 5e-3)
        #expect(abs(segs[1].hi - 150.4075) < 5e-3)
        #expect(abs((segs[1].lo - segs[0].hi) - 13.8087) < 5e-3, "第一禁带")
        // P 减小 → 禁带变窄（布拉格反射减弱）
        let (_, allowed2, _) = scan(P: 0.2)
        let segs2 = KronigPenneyMath.segments(allowed: allowed2, E: E)
        #expect(segs2.count >= 2 && segs2[1].lo - segs2[0].hi < segs[1].lo - segs[0].hi,
                "P=0.2 禁带应窄于 P=1")
    }

    @Test("gradient：内部中心差分、两端单侧（np.gradient edge_order=1 语义）")
    func gradientKernel() {
        #expect(KronigPenneyMath.gradient([1, 2, 4, 8, 16]) == [1, 1.5, 3, 6, 8])
        #expect(KronigPenneyMath.gradient([3]) == [3])
        #expect(KronigPenneyMath.gradient([5, 7]) == [2, 2])
    }

    @Test("percentile：线性插值法（np.nanpercentile 默认 linear）")
    func percentileKernel() {
        #expect(KronigPenneyMath.percentile([1, 2, 3, 4], 0) == 1)
        #expect(KronigPenneyMath.percentile([1, 2, 3, 4], 50) == 2.5)
        #expect(KronigPenneyMath.percentile([1, 2, 3, 4], 100) == 4)
        #expect(KronigPenneyMath.percentile([1, 2, 3, 4], 25) == 1.75)
        #expect(KronigPenneyMath.percentile([1, 2, 3, 4], 98) == 3.94)
    }

    @Test("fixture 曲线对拍：g(E)（512 点，完整重算序列最近邻匹配，cap = 26.5794）")
    func dosFixtureMatch() throws {
        let (E, allowed, kRed) = scan()
        let segs = KronigPenneyMath.segments(allowed: allowed, E: E)
        var eAll: [Double] = [], gAll: [Double] = []
        for seg in segs {
            var es: [Double] = [], ks: [Double] = []
            for i in E.indices where E[i] >= seg.lo - 1e-9 && E[i] <= seg.hi + 1e-9 {
                es.append(E[i]); ks.append(kRed[i])
            }
            guard es.count >= 4 else { continue }
            let g = KronigPenneyMath.dosInSegment(E: es, k: ks)
            eAll.append(contentsOf: es)
            gAll.append(contentsOf: g)
        }
        var finite: [(e: Double, g: Double)] = []
        for i in gAll.indices where gAll[i].isFinite { finite.append((eAll[i], gAll[i])) }
        #expect(finite.count == 16263, "有限值点数 = \(finite.count)（Python 端 16263）")
        let cap = KronigPenneyMath.percentile(finite.map(\.g).sorted(), 98)
        #expect(abs(cap - 26.579424) < 1e-5, "cap98 = \(cap)")
        let gClip = finite.map { min($0.g, cap) }

        let fx = try ModuleFixture.load("晶体中的电子_能带与KronigPenney__v2022")
        let (x, y) = try fx.line(1, 0, index: 0)
        #expect(x.count == 512)
        #expect(abs(x.first! - 6.50442) < 1e-3 && abs(x.last! - 150.408) < 1e-3)
        #expect(abs(y.max()! - cap) < 1e-5, "截断峰 = \(y.max()!) vs cap = \(cap)")
        #expect(y.min()! >= 0, "g ≥ 0")
        // eAll 单调递增 → 双指针最近邻匹配（网格间距 8e-3 eV，x 舍入 ~1e-7 无歧义）
        var j = 0
        for i in x.indices {
            while j + 1 < finite.count,
                  abs(finite[j + 1].e - x[i]) < abs(finite[j].e - x[i]) { j += 1 }
            #expect(abs(finite[j].e - x[i]) < 1e-5, "第 \(i) 点 x 漂移：\(finite[j].e) vs \(x[i])")
            #expect(abs(gClip[j] - y[i]) < 1e-6 * max(1, abs(y[i])),
                    "g(E) 第 \(i) 点：\(gClip[j]) vs \(y[i])")
        }
    }

    @Test("compute 输出结构：2 图（能带 V 形折线 + DOS）+ 摘要 4 条 + 理论卡")
    func computeStructure() async throws {
        let module = KronigPenneyModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        guard case .lineSeries(let c0) = result.charts[0],
              case .lineSeries(let c1) = result.charts[1] else {
            Issue.record("应为 2 lineSeries"); return
        }
        // 默认 P=1：恰 2 条带系列，±k 双支，k ∈ [−π/a, π/a]
        #expect(c0.series.count == 2)
        for s in c0.series {
            #expect(s.points.count > 100)
            for p in s.points { #expect(abs(p.x) <= .pi + 1e-9) }
        }
        #expect(c0.referenceLines.count == 2)
        #expect(c1.series[0].points.count >= 500, "DOS 降采样 ≥ 500 点")
        #expect(result.summary.count == 4)
        #expect(result.summary[1].value.contains("13.80"), "第一禁带 = \(result.summary[1].value)")
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("秒级档预算：compute 中位 < 2000 ms（20000 点扫描 + DOS）")
    func computeBudget() async throws {
        let module = KronigPenneyModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 2000)
    }
}
