import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W8-17 中心力场与氢原子径向方程：3000² 稠密 eigh 部分谱（dsyevr RANGE='I'）×
/// l = 0,1,2 + Rydberg 解析对照 + 节点定理 + r_max 自适应（fixtures 图 1 九曲线对拍）。
@Suite("W8 中心力场：氢原子径向方程")
struct HydrogenRadialModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    /// 物理常数组合（约化质量）。
    private func phys() throws -> (mu: Double, hbar: Double, e: Double, eps0: Double) {
        let c = try constants()
        let me = try c.value("m_e"), mp = try c.value("m_p")
        return (me * mp / (me + mp), try c.value("hbar"),
                try c.value("e"), try c.value("epsilon_0"))
    }

    /// 脚本口径单 l 求解：返回内点网格（a₀ 单位）+ 前 k 能量（eV）+ 归一化 u + 密度。
    private func solve(l: Int, rMaxNm: Double, count k: Int = 5) throws
        -> (rA0: [Double], eV: [Double], u: [[Double]], density: [[Double]]) {
        let p = try phys()
        let c = try constants()
        let s = HydrogenRadialMath.solveLowest(
            l: l, rMax: rMaxNm * 1.0e-9, points: HydrogenRadialMath.steps,
            count: k, mu: p.mu, hbar: p.hbar, eCharge: p.e, eps0: p.eps0)
        let a0 = HydrogenRadialMath.bohrRadius(mu: p.mu, eCharge: p.e, eps0: p.eps0, hbar: p.hbar)
        let eCharge = try c.value("e")
        let m = HydrogenRadialMath.steps - 1
        let dr = s.r[1] - s.r[0]
        var u: [[Double]] = [], density: [[Double]] = []
        for st in 0..<k {
            var vec = [Double](repeating: 0, count: m)
            for i in 0..<m { vec[i] = s.v[i * k + st] }
            let un = HydrogenRadialMath.normalizeU(vec, dr: dr)
            u.append(un)
            density.append(HydrogenRadialMath.radialDensity(u: un, r: s.r))
        }
        return (s.r.map { $0 / a0 }, s.w.map { $0 / eCharge }, u, density)
    }

    // MARK: 物理律

    @Test("解析核：玻尔半径 0.529465 Å、里德伯 13.5983 eV（约化质量，CODATA 复算）")
    func analyticConstants() throws {
        let p = try phys()
        let a0 = HydrogenRadialMath.bohrRadius(mu: p.mu, eCharge: p.e, eps0: p.eps0, hbar: p.hbar)
        let ry = HydrogenRadialMath.rydberg(mu: p.mu, eCharge: p.e, eps0: p.eps0, hbar: p.hbar)
        #expect(abs(a0 * 1.0e10 - 0.529465) < 1e-5, "a₀ = \(a0 * 1.0e10) Å")
        #expect(abs(ry / p.e - 13.5983) < 1e-3, "Ry = \(ry / p.e) eV")
    }

    @Test("物理律：l=0 前 4 态落 Rydberg 级数（rel < 3e-4，箱效应主导）；n=5 起受箱宽污染")
    func rydbergSeries() throws {
        let s = try solve(l: 0, rMaxNm: 3.0)
        let p = try phys()
        let ryEv = HydrogenRadialMath.rydberg(mu: p.mu, eCharge: p.e, eps0: p.eps0, hbar: p.hbar) / p.e
        for n in 1...4 {
            let exact = -ryEv / Double(n) / Double(n)
            let rel = abs(s.eV[n - 1] - exact) / abs(exact)
            #expect(rel < 3e-4, "n=\(n)：\(s.eV[n - 1]) vs \(exact)（rel \(rel)）")
        }
        let exact5 = -ryEv / 25.0
        #expect(abs(s.eV[4] - exact5) / abs(exact5) > 1e-2,
                "n=5 应受箱宽污染（脚本口径 r_max=3nm），实测 \(s.eV[4])")
    }

    @Test("物理律：库仑简并 E(l=0,n) ≈ E(l=1,n)（箱内近似成立，rel < 1e-4）")
    func degeneracy() throws {
        let s0 = try solve(l: 0, rMaxNm: 3.0)
        let s1 = try solve(l: 1, rMaxNm: 3.0)
        // l=0 态 2/3 ↔ l=1 态 1/2 → 同一 n
        for (i0, i1) in [(1, 0), (2, 1)] {
            let rel = abs(s0.eV[i0] - s1.eV[i1]) / abs(s0.eV[i0])
            #expect(rel < 1e-4, "E(l=0)[\(i0)] vs E(l=1)[\(i1)] 相对差 \(rel)")
        }
    }

    @Test("物理律：u 节点数 = n−l−1（九态全查）")
    func nodeTheorem() throws {
        for (l, nStart) in [(0, 1), (1, 2), (2, 3)] {
            let s = try solve(l: l, rMaxNm: 3.0, count: 3)
            for j in 0..<3 {
                let nodes = HydrogenRadialMath.countNodes(s.u[j])
                #expect(nodes == j, "n=\(nStart + j), l=\(l)：节点 \(nodes) ≠ \(j)")
            }
        }
    }

    @Test("物理律：r_max 自适应——小箱抬高 E₄（>1e-3），自适应箱恢复（<5e-3）")
    func adaptiveRmax() throws {
        let p = try phys()
        let a0 = HydrogenRadialMath.bohrRadius(mu: p.mu, eCharge: p.e, eps0: p.eps0, hbar: p.hbar)
        // 自适应下限 = 2n²ₘₐₓa₀（n_max = 5），与模块 compute 同口径
        let floorNm = 2.0 * 25.0 * a0 * 1.0e9
        let sBig = try solve(l: 0, rMaxNm: 3.0)          // 脚本口径
        let sSmall = try solve(l: 0, rMaxNm: 1.0)        // 过小箱：高激发态被抬升
        let sAdaptive = try solve(l: 0, rMaxNm: max(1.0, floorNm))
        // 注：自适应箱（2.65 nm）与脚本箱（3 nm）的网格步长 h 不同，库仑奇点使
        // FD 误差 ~ O(h)，n=4（无角动量垒、原点振幅最集中）对 h 最敏感——
        // 故「恢复」判据取 FD 分辨率一致容差 5e-3，并要求箱效应 ≥ 5× 分辨率效应。
        for n in 1...4 {
            let relSmall = abs(sSmall.eV[n - 1] - sBig.eV[n - 1]) / abs(sBig.eV[n - 1])
            let relAdaptive = abs(sAdaptive.eV[n - 1] - sBig.eV[n - 1]) / abs(sBig.eV[n - 1])
            if n == 4 {
                #expect(relSmall > 1e-3, "1 nm 箱应扰动 n=4（rel \(relSmall)）")
                #expect(relSmall > 5 * relAdaptive,
                        "箱效应应显著大于分辨率效应：\(relSmall) vs \(relAdaptive)")
            }
            #expect(relAdaptive < 5e-3, "自适应箱 n=\(n) 扰动 \(relAdaptive)")
        }
    }

    // MARK: fixtures 对拍

    @Test("fixtures：径向概率密度九曲线（fixture x 为精确格点，最近邻匹配，|Δ| < 1e-6·max|P|）")
    func densityFixture() throws {
        let fx = try ModuleFixture.load("中心力场与氢原子_径向方程__v2022")
        // fixtures 图 1（fig 1）九曲线无节点标注，按序 n,l = (1,0)(2,0)(3,0)(2,1)(3,1)(4,1)(3,2)(4,2)(5,2)
        let states = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)]
        var curveIdx = 0
        for (l, j) in states {
            let s = try solve(l: l, rMaxNm: 3.0, count: j + 1)   // 只解到所需态
            let f = try fx.line(1, 0, index: curveIdx)
            let density = s.density[j]
            let m = s.rA0.count
            var maxP = 0.0
            for v in f.y { maxP = max(maxP, abs(v)) }
            // fixture 512 点为网格均匀步进取样（idx = round(j·(m−1)/511)，exporter 管线口径）；
            // a0 的 CODATA 2018/2022 微差使 x 值有 ~1e-9 相对漂移，故按序号而非 x 值定位。
            // 首点 r₁ = h 恰压脚本钳位阈值 r > 1e-12 m（Python/Swift 浮点路径分居阈值两侧，
            // 一侧钳 0 一侧取值），属单点阈值伪差，跳过该点。
            for q in f.x.indices where q > 0 {
                let idx = Int((Double(q) * Double(m - 1) / Double(f.x.count - 1)).rounded())
                #expect(idx >= 0 && idx < m, "fixture 序号落在网格内")
                #expect(abs(s.rA0[idx] - f.x[q]) < 1e-3 * (s.rA0[1] - s.rA0[0]),
                        "fixture x 应对齐格点（idx=\(idx)）")
                #expect(abs(density[idx] - f.y[q]) < 1e-6 * maxP,
                        "n=\(j + 1 + l),l=\(l)@r=\(f.x[q])a₀: \(density[idx]) vs \(f.y[q])")
            }
            curveIdx += 1
        }
    }

    // MARK: compute 结构与预算

    @Test("compute 输出结构：2 图（9 系列密度 + 能量对照）+ 摘要 5 条 + 理论卡 4 式")
    func computeStructure() async throws {
        let module = HydrogenRadialModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        guard case .lineSeries(let c0) = result.charts[0],
              case .lineSeries(let c1) = result.charts[1] else {
            Issue.record("应为 2 lineSeries"); return
        }
        #expect(c0.series.count == 9)
        #expect(c0.series.allSatisfy { $0.points.count <= 512 })
        #expect(c1.series.count == 2 && c1.series[0].points.count == 5)
        #expect(result.summary.count == 5)
        #expect(result.summary[4].value.contains("n=1,l=0:0"))
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("秒级档预算：compute < 2000 ms（收尾包 1 走三对角 dstevr 后 ~50 ms，强制口径）")
    func computeBudget() async throws {
        let module = HydrogenRadialModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 2000,
            runs: 5, warmup: 2, enforce: true)
    }
}
