import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W8-13 薛定谔方程有限差分势阱：800² 三对角 eigh 部分谱 + E_n ∝ n² 解析对照 +
/// 节点定理（Python 端无 [check] 断言，fixtures 图曲线/能量点复算对拍）。
@Suite("W8 薛定谔方程：有限差分势阱")
struct SchrodingerFDWellModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    /// 脚本口径求解（L nm → m，N = 800，前 6 对），返回 eV 能量 + 网格归一化 ψ。
    private func solve(widthNm L: Double) throws -> (eV: [Double], analyticEv: [Double],
                                                     psi: [[Double]], xNm: [Double]) {
        let c = try constants()
        let hbar = try c.value("hbar"), me = try c.value("m_e"), e = try c.value("e")
        let w = L * 1.0e-9
        let n = SchrodingerFDWellMath.gridPoints
        let dx = w / Double(n + 1)
        let s = SchrodingerFDWellMath.solveLowest(
            widthMeters: w, points: n, count: 6, hbar: hbar, mass: me)
        let eV = s.w.map { $0 / e }
        let analytic = SchrodingerFDWellMath.analyticEnergies(
            widthMeters: w, count: 6, hbar: hbar, mass: me).map { $0 / e }
        var psi: [[Double]] = []
        for st in 0..<6 {
            var vec = [Double](repeating: 0, count: n)
            for i in 0..<n { vec[i] = s.v[i * 6 + st] }
            psi.append(SchrodingerFDWellMath.normalizeOnGrid(vec, dx: dx))
        }
        let xNm = (1...n).map { Double($0) * dx * 1.0e9 }
        return (eV, analytic, psi, xNm)
    }

    // MARK: 物理律

    @Test("解析核：L = 1 nm 基态 = π²ħ²/2mL² = 0.376030 eV（CODATA 精确复算）")
    func analyticKernel() throws {
        let c = try constants()
        let hbar = try c.value("hbar"), me = try c.value("m_e"), eC = try c.value("e")
        let e1 = SchrodingerFDWellMath.analyticEnergies(
            widthMeters: 1.0e-9, count: 1, hbar: hbar, mass: me)[0] / eC
        #expect(abs(e1 - 0.376030) < 1e-6, "E₁ 解析 = \(e1)")
        // E ∝ n² 严格成立
        let all = SchrodingerFDWellMath.analyticEnergies(
            widthMeters: 1.0e-9, count: 6, hbar: hbar, mass: me).map { $0 / eC }
        for n in 1...6 {
            #expect(abs(all[n - 1] / all[0] - Double(n) * Double(n)) < 1e-12)
        }
    }

    @Test("物理律：FD 数值谱从下方逼近解析（n≤6 rel < 5e-5）且节点数 = n−1")
    func spectralConvergenceAndNodes() throws {
        let (eV, analytic, psi, _) = try solve(widthNm: 1.0)
        for i in 0..<6 {
            let rel = abs(eV[i] - analytic[i]) / analytic[i]
            #expect(rel < 5e-5, "E_\(i + 1) 相对误差 \(rel)")
            #expect(eV[i] < analytic[i], "FD 刚度误差使数值谱从下方逼近（sin 展开二阶项）")
            let nodes = SchrodingerFDWellMath.countNodes(psi[i])
            #expect(nodes == i, "n=\(i + 1) 节点数 \(nodes) ≠ \(i)")
        }
        // E_n ∝ n²（FD 层面，1e-3 容差）
        for n in 2...6 {
            #expect(abs(eV[n - 1] / eV[0] - Double(n) * Double(n)) < 5e-3)
        }
    }

    @Test("阱宽标度：E₁(L) ∝ 1/L²（L 翻倍 → E₁/4）")
    func widthScaling() throws {
        let (e1, _, _, _) = try solve(widthNm: 1.0)
        let (e2, _, _, _) = try solve(widthNm: 2.0)
        #expect(abs(e1[0] / e2[0] - 4.0) < 1e-3,
                "E₁(1nm)/E₁(2nm) = \(e1[0] / e2[0])")
    }

    // MARK: fixtures 对拍

    @Test("fixtures：能量点（解析 + 数值，6 点，rel < 1e-8）")
    func energyFixture() throws {
        let (eV, analytic, _, _) = try solve(widthNm: 1.0)
        let fx = try ModuleFixture.load("薛定谔方程_有限差分势阱__v2022")
        let analyticFix = try fx.line(0, 1, index: 0)   // "analytic n^2 pi^2 ..."
        let numericFix = try fx.line(0, 1, index: 1)    // "finite-difference numeric"
        #expect(analyticFix.x.count == 6 && numericFix.x.count == 6)
        for i in 0..<6 {
            #expect(abs(analyticFix.x[i] - Double(i + 1)) < 1e-12, "x = n")
            expectPointwiseClose([analytic[i]], [analyticFix.y[i]], tolerance: 1e-8,
                                 "解析 E_\(i + 1)")
            expectPointwiseClose([eV[i]], [numericFix.y[i]], tolerance: 1e-8,
                                 "数值 E_\(i + 1)")
        }
    }

    @Test("fixtures：本征函数 ψ₁–ψ₄（fixture x 为精确格点，最近邻匹配，|Δ| < 1e-6·max|ψ|）")
    func eigenfunctionFixture() throws {
        let (_, _, psi, xNm) = try solve(widthNm: 1.0)
        let dxNm = xNm[1] - xNm[0]
        let fx = try ModuleFixture.load("薛定谔方程_有限差分势阱__v2022")
        for s in 0..<4 {
            let f = try fx.line(0, 0, index: s)   // "n=\(s+1)"（含 1.2s 偏置）
            let offset = 1.2 * Double(s)
            let fy = f.y.map { $0 - offset }      // 还原 ψ
            // 符号对齐（特征向量符号约定因 LAPACK 实现/版本而异）
            var dot = 0.0
            for j in f.x.indices {
                let idx = Int((f.x[j] / dxNm).rounded()) - 1
                dot += psi[s][idx] * fy[j]
            }
            let sg: Double = dot < 0 ? -1 : 1
            var maxPsi = 0.0
            for v in fy { maxPsi = max(maxPsi, abs(v)) }
            for j in f.x.indices {
                let idx = Int((f.x[j] / dxNm).rounded()) - 1
                #expect(abs(Double(idx + 1) * dxNm - f.x[j]) < 1e-9,
                        "fixture x 应为精确格点（idx=\(idx)）")
                #expect(abs(sg * psi[s][idx] - fy[j]) < 1e-6 * maxPsi,
                        "ψ_\(s + 1)@x=\(f.x[j]): \(sg * psi[s][idx]) vs \(fy[j])")
            }
        }
    }

    // MARK: compute 结构与预算

    @Test("compute 输出结构：2 图（4+2 系列）+ 摘要 4 条 + 理论卡 4 式")
    func computeStructure() async throws {
        let module = SchrodingerFDWellModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        guard case .lineSeries(let c0) = result.charts[0],
              case .lineSeries(let c1) = result.charts[1] else {
            Issue.record("应为 2 lineSeries"); return
        }
        #expect(c0.series.count == 4)
        #expect(c0.series.allSatisfy { $0.points.count == 800 })
        #expect(c1.series.count == 2 && c1.series[0].points.count == 6)
        #expect(result.summary.count == 4)
        #expect(result.summary[2].value == "0, 1, 2, 3", "节点计数 = \(result.summary[2].value)")
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("秒级档预算：compute < 2000 ms（800² 部分谱，强制口径）")
    func computeBudget() async throws {
        let module = SchrodingerFDWellModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 2000)
    }
}
