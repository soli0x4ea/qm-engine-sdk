import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W8-14 一维势场：矩形势垒精确透射 + 有限深势阱超越方程二分求根。
/// 脚本无 [check] 断言；T(E)/T(a) 闭式曲线逐点对拍，波函数曲线对拍
/// （只对拍真根态：脚本扫描把 tan 极点跳变误判为根——**基准脚本疑点**，见周报）。
@Suite("W8 一维势场：方势垒与有限深势阱")
struct BarrierWellModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    private func phys() throws -> (hbar: Double, me: Double, e: Double) {
        let c = try constants()
        return (try c.value("hbar"), try c.value("m_e"), try c.value("e"))
    }

    // MARK: 势垒（物理律）

    @Test("T(E) 核：脚本锚点 T(0.3 eV, a=1 nm) = 6.3555e-4；厚垒极限 6.3572e-4")
    func transmissionAnchor() throws {
        let p = try phys()
        let v0 = 1.0 * p.e
        let t = BarrierWellMath.transmission(E: 0.3 * p.e, v0: v0, width: 1.0e-9,
                                             mass: p.me, hbar: p.hbar)
        let tl = BarrierWellMath.thickLimit(E: 0.3 * p.e, v0: v0, width: 1.0e-9,
                                            mass: p.me, hbar: p.hbar)
        #expect(abs(t - 6.3555e-4) / 6.3555e-4 < 1e-4, "T = \(t)")
        #expect(abs(tl - 6.3572e-4) / 6.3572e-4 < 1e-4, "厚垒极限 = \(tl)")
        // 越域 NaN（脚本口径）
        #expect(BarrierWellMath.transmission(E: 0.0, v0: v0, width: 1e-9,
                                             mass: p.me, hbar: p.hbar).isNaN)
        #expect(BarrierWellMath.transmission(E: v0, v0: v0, width: 1e-9,
                                             mass: p.me, hbar: p.hbar).isNaN)
    }

    @Test("物理律：E > V₀ 越垒共振 T = 1 落在 kL = nπ（n=3 精确复核）")
    func transmissionResonance() throws {
        let p = try phys()
        let v0 = 1.0 * p.e
        let a = 1.0e-9
        // kL = 3π → E = ħ²(3π/a)²/2m
        let k = 3.0 * .pi / a
        let eRes = p.hbar * p.hbar * k * k / (2.0 * p.me)
        #expect(eRes > v0, "共振能须在垒上")
        let t = BarrierWellMath.transmissionAbove(E: eRes, v0: v0, width: a,
                                                  mass: p.me, hbar: p.hbar)
        #expect(abs(t - 1.0) < 1e-9, "共振 T = \(t)")
        // 反共振点（kL = 3.5π）：T 明显 < 1（薄垒共振宽，0.99 量级）
        let kAnti = 3.5 * .pi / a
        let eAnti = p.hbar * p.hbar * kAnti * kAnti / (2.0 * p.me)
        let tAnti = BarrierWellMath.transmissionAbove(E: eAnti, v0: v0, width: a,
                                                      mass: p.me, hbar: p.hbar)
        #expect(tAnti < 0.999, "反共振 T = \(tAnti)")
        #expect(tAnti < t - 1e-6)
    }

    // MARK: 有限深势阱（物理律 + 脚本疑点）

    @Test("物理律：束缚态个数 = ⌈2z₀/π⌉（V₀=0.5 eV/L=2 nm → 3；V₀=1.5 eV → 4；单调）")
    func boundStateCount() throws {
        let p = try phys()
        let s1 = BarrierWellMath.boundStates(v0: 0.5 * p.e, width: 2.0e-9,
                                             mass: p.me, hbar: p.hbar)
        let s2 = BarrierWellMath.boundStates(v0: 1.5 * p.e, width: 2.0e-9,
                                             mass: p.me, hbar: p.hbar)
        #expect(s1.count == BarrierWellMath.analyticBoundCount(
            v0: 0.5 * p.e, width: 2.0e-9, mass: p.me, hbar: p.hbar))
        #expect(s1.count == 3, "实测 \(s1.count)")
        #expect(s2.count == 4 && s2.count > s1.count, "实测 \(s2.count)")
        // 宇称按能量序交替：even/odd/even
        #expect(s1[0].even && !s1[1].even && s1[2].even)
    }

    @Test("物理律：束缚态能量（V₀=0.5 eV, L=2 nm）= 脚本真值三态（rel < 1e-4）")
    func boundStateEnergies() throws {
        let p = try phys()
        let s = BarrierWellMath.boundStates(v0: 0.5 * p.e, width: 2.0e-9,
                                            mass: p.me, hbar: p.hbar)
        let truthEv = [-0.44277, -0.27821, -0.04630]   // Python 复算（剔除伪根后）
        for (i, eTrue) in truthEv.enumerated() {
            let eEv = s[i].E / p.e
            #expect(abs(eEv - eTrue) / abs(eTrue) < 1e-4, "E_\(i + 1) = \(eEv) vs \(eTrue)")
        }
        #expect(s[0].E > -0.5 * p.e && s[2].E < 0, "能量 ∈ (−V₀, 0)")
    }

    @Test("基准脚本疑点：原始扫描把 tan 极点误判为根（5 个），剔除后剩 3 个真根")
    func spuriousRoots() throws {
        let p = try phys()
        let v0 = 0.5 * p.e, l = 2.0e-9
        let kmax = sqrt(2.0 * p.me * v0) / p.hbar
        let raw = BarrierWellMath.rawRoots(v0: v0, width: l, mass: p.me, hbar: p.hbar)
        let valid = raw.filter { BarrierWellMath.isValidRoot($0.k, even: $0.even, kmax: kmax, width: l) }
        #expect(raw.count == 5, "脚本口径原始根数 \(raw.count)")
        #expect(valid.count == 3, "真根数 \(valid.count)")
        // 伪根坐落 tan 极点 kL/2 = π/2, π
        let spurious = raw.filter { !BarrierWellMath.isValidRoot($0.k, even: $0.even, kmax: kmax, width: l) }
        for s in spurious {
            let phase = s.k * l / 2.0
            let nearPole = abs(phase - .pi / 2) < 1e-4 || abs(phase - .pi) < 1e-4
            #expect(nearPole, "伪根相位 \(phase) 应在 tan/cot 极点")
        }
    }

    // MARK: fixtures 对拍

    @Test("fixtures：T(E) 400 点逐点（rel < 1e-8）+ T(a)/厚垒极限 400 点逐点")
    func barrierFixture() throws {
        let p = try phys()
        let v0 = 1.0 * p.e
        let fx = try ModuleFixture.load("一维势场_方势垒与有限深势阱__v2022")
        // 图 0 轴 0：T(E)，E ∈ [0.02, 0.98]·V₀
        let (ex, ey) = try fx.line(0, 0, index: 0)
        #expect(ex.count == 400)
        for i in ex.indices {
            let t = BarrierWellMath.transmission(E: ex[i] * v0, v0: v0, width: 1.0e-9,
                                                 mass: p.me, hbar: p.hbar)
            expectPointwiseClose([t], [ey[i]], tolerance: 1e-8, "T(E=\(ex[i])V₀)")
        }
        // 图 0 轴 1：T(a) 精确 + 厚垒极限（E = 0.3 eV）
        let (ax, ay) = try fx.line(0, 1, index: 0)
        let (lx, ly) = try fx.line(0, 1, index: 1)
        // T(a) 含 e^(−2κa) 陡峭指数：fixture x 9 位存储的 δa/a≈5e-9 经 2κa（≤23）放大，
        // 存储精度封顶 ~2e-7（同 W7 洛伦兹斜率放大口径），故此二曲线容差放宽到 2e-7
        for i in ax.indices {
            let am = ax[i] * 1.0e-9
            let t = BarrierWellMath.transmission(E: 0.3 * p.e, v0: v0, width: am,
                                                 mass: p.me, hbar: p.hbar)
            let tl = BarrierWellMath.thickLimit(E: 0.3 * p.e, v0: v0, width: am,
                                                mass: p.me, hbar: p.hbar)
            expectPointwiseClose([t], [ay[i]], tolerance: 2e-7, "T(a=\(ax[i])nm)")
            expectPointwiseClose([tl], [ly[i]], tolerance: 2e-7, "厚垒极限(a=\(ax[i])nm)")
        }
        #expect(lx.count == ax.count)
    }

    @Test("fixtures：真根态波函数（fixture n=1/n=5 为真态，符号对齐 |Δ| < 2e-6·max|ψ|）")
    func wavefunctionFixture() throws {
        let p = try phys()
        let v0 = 0.5 * p.e, l = 2.0e-9
        let states = BarrierWellMath.boundStates(v0: v0, width: l, mass: p.me, hbar: p.hbar)
        let fx = try ModuleFixture.load("一维势场_方势垒与有限深势阱__v2022")
        // fixture 图 1 轴 0：n=1（idx0 → 真态 0）、n=5（idx4 → 真态 2，真偶宇称）
        // 容差 2e-6：fixture x 是 2000 点网格索引重映射的抖动网格，梯形法归一的
        // 剩余网格误差包络 ~7e-7（n=3 峰曲率更大），取 3 倍裕度。
        for (fixIdx, stateIdx) in [(0, 0), (4, 2)] {
            let f = try fx.line(1, 0, index: fixIdx)
            let st = states[stateIdx]
            let offset = 0.6 * Double(fixIdx)
            let fy = f.y.map { $0 - offset }
            let xsM = f.x.map { $0 * 1.0e-9 }
            let psi = BarrierWellMath.wavefunction(k: st.k, v0: v0, width: l, even: st.even,
                                                   xs: xsM, mass: p.me, hbar: p.hbar)
            var dot = 0.0
            for i in xsM.indices { dot += psi[i] * fy[i] }
            let sg: Double = dot < 0 ? -1 : 1
            var maxPsi = 0.0
            for v in fy { maxPsi = max(maxPsi, abs(v)) }
            for i in xsM.indices {
                #expect(abs(sg * psi[i] - fy[i]) < 2e-6 * maxPsi,
                        "n=\(stateIdx + 1)@x=\(f.x[i]): \(sg * psi[i]) vs \(fy[i])")
            }
        }
    }

    // MARK: compute 结构与预算

    @Test("compute 输出结构：4 图 + 摘要 4 条 + 理论卡 4 式（默认参数 3 束缚态）")
    func computeStructure() async throws {
        let module = BarrierWellModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 4)
        guard case .lineSeries(let c0) = result.charts[0],
              case .lineSeries(let c1) = result.charts[1],
              case .lineSeries(let c2) = result.charts[2],
              case .lineSeries(let c3) = result.charts[3] else {
            Issue.record("应为 4 lineSeries"); return
        }
        #expect(c0.series[0].points.count == 400)
        #expect(c1.series.count == 2 && c1.series[0].points.count == 400)
        #expect(c2.series.count == 3, "默认阱深 3 个真束缚态，实测 \(c2.series.count)")
        #expect(c2.series[0].points.count == 2000)
        #expect(c3.series.count == 3 && c3.referenceLines.count == 2)
        #expect(result.summary.count == 4)
        #expect(result.summary[1].value.contains("3（解析 ⌈2z₀/π⌉ = 3）"))
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("秒级档预算：compute < 2000 ms（200000 点扫描 ×2 + 二分，强制口径）")
    func computeBudget() async throws {
        let module = BarrierWellModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 2000)
    }
}
