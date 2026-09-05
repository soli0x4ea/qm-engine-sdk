import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W9-38 狄拉克方程：能谱。
/// 4×4 实对称 H(p) 400 点扫描（dsyevr）+ 解析包络 fixture 逐点对拍
/// + 谱对称/二重简并/质量隙/非相对论极限物理律。
@Suite("W9 狄拉克方程：能谱")
struct DiracSpectrumModuleTests {

    private func constants() throws -> ConstantsSet {
        let c = try ConstantsSet.load(.v2022)
        return c
    }

    private func meC2() throws -> Double { try constants().value("me_c2_MeV") }

    // MARK: 计算核（物理律）

    @Test("脚本锚点：p=0 → ±mc²（各 2 重）= ±0.510999；p=1 → ±1.122996")
    func scriptAnchors() throws {
        let m = try meC2()
        let e0 = AlgebraCore.eigh(DiracSpectrumMath.hamiltonian(0.0, meC2: m), n: 4).w
        #expect(abs(e0[3] - m) < 1e-9 && abs(e0[2] - m) < 1e-9, "正能 2 重：\(e0[2]), \(e0[3])")
        #expect(abs(e0[1] + m) < 1e-9 && abs(e0[0] + m) < 1e-9, "负能 2 重：\(e0[0]), \(e0[1])")
        #expect(abs(m - 0.510999) < 5e-6, "mc² = \(m)（脚本打印 0.510999）")

        let e1 = AlgebraCore.eigh(DiracSpectrumMath.hamiltonian(1.0, meC2: m), n: 4).w
        #expect(abs(e1[3] - 1.122996) < 5e-6, "E₊(p=1) = \(e1[3])")
        #expect(abs(e1[0] + 1.122996) < 5e-6, "E₋(p=1) = \(e1[0])")
    }

    @Test("物理律：谱对称 E↔−E 与解析 |E| rel < 1e-12（400 点全扫描）")
    func spectralSymmetryAndAnalytic() throws {
        let m = try meC2()
        let pGrid = Num.linspace(0, 3.0, count: 400)
        let eNum = DiracSpectrumMath.scan(pGrid: pGrid, meC2: m)
        for (i, p) in pGrid.enumerated() {
            let w = eNum[i]
            // 谱对称：升序 w，w[3] = −w[0]、w[2] = −w[1]
            #expect(abs(w[3] + w[0]) < 1e-12 * max(1.0, abs(w[3])))
            #expect(abs(w[2] + w[1]) < 1e-12 * max(1.0, abs(w[2])))
            // 解析：正能支与 |E| = √(p²+m²) 一致（两支简并）
            let analytic = DiracSpectrumMath.analyticEnergy(p, meC2: m)
            let rel = abs(w[3] - analytic) / analytic
            #expect(rel < 1e-12, "p = \(p): E = \(w[3]) vs \(analytic)")
        }
    }

    @Test("物理律：非相对论极限 E₊ − mc² ≈ p²/2m（p = 0.01，rel < 5e-4）；质量隙 1.022 MeV")
    func nonrelativisticLimitAndGap() throws {
        let m = try meC2()
        let p = 0.01
        let e = AlgebraCore.eigh(DiracSpectrumMath.hamiltonian(p, meC2: m), n: 4).w
        let kinetic = e[3] - m
        let classical = p * p / (2.0 * m)
        // 下一阶修正 −p⁴/(8m³c⁴) 的相对量级 ≈ p²/(4m²c⁴) ≈ 1e-4（p=0.01 MeV/c）
        #expect(abs(kinetic - classical) / classical < 5e-4,
                "动能 \(kinetic) vs p²/2m \(classical)")
        #expect(abs(DiracSpectrumMath.massGap(meC2: m) - 1.022) < 5e-4)
    }

    // MARK: fixture 逐点对拍（解析包络 400 点原生网格）

    @Test("fixture 对拍：解析包络 |E| 400 点（rel < 1e-8）")
    func analyticEnvelopeFixture() throws {
        let fixture = try ModuleFixture.load("狄拉克方程_能谱__v2022")
        let fx = try fixture.line(0, 0, label: "analytic")
        let m = try meC2()
        let ours = fx.x.map { DiracSpectrumMath.analyticEnergy($0, meC2: m) }
        expectPointwiseClose(ours, fx.y, tolerance: 1e-8, "解析包络")
    }

    // MARK: 模块契约

    @Test("模块契约：单图 6 序列各 400 点 + 摘要 5 卡")
    func moduleContract() async throws {
        let module = DiracSpectrumModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())

        #expect(result.charts.count == 1)
        guard case .lineSeries(let chart) = result.charts[0] else {
            Issue.record("图 1 应为曲线"); return
        }
        #expect(chart.series.count == 6)
        #expect(chart.series.allSatisfy { $0.points.count == 400 })

        #expect(result.summary.count == 5)
        #expect(result.summary[2].value.contains("2mc²"))
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("秒级档预算：compute < 2000 ms（400 × 4×4 dsyevr）")
    func computeBudget() async throws {
        let module = DiracSpectrumModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 2000, runs: 5, warmup: 2)
    }
}
