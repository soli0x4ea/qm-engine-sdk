import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W8-19 微扰论与变分法：40² 谐振子基底二阶微扰 vs 精确对角化 + 方势阱斜坡
/// （一阶高估 ~2×）+ 氦单参变分（上界定理）。脚本无 [check] 断言，
/// 数值锚点取 Python 复算打印值；氦 E(α) 扫描曲线 fixtures 对拍。
@Suite("W8 微扰论与变分法")
struct PerturbationVariationModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    // MARK: (A1) 非谐振子

    @Test("基底核：⟨X⁴⟩₀ = 3/4（谐振子基态四阶矩恒等）")
    func x4MatrixElement() {
        let (e0, x4) = PerturbationVariationMath.anharmonicBasis(N: 40)
        #expect(e0[0] == 0.5 && e0[1] == 1.5)
        #expect(abs(x4[0] - 0.75) < 1e-12, "⟨X⁴⟩₀ = \(x4[0])")
        // X⁴ 对称
        #expect(abs(x4[0 * 40 + 3] - x4[3 * 40 + 0]) < 1e-12)
    }

    @Test("物理律：二阶微扰 (0+1+2) 与精确对角化一致到 λ⁴ 量级（λ=0.05/0.10 脚本锚点）")
    func secondOrderAnchors() {
        let (e0, x4) = PerturbationVariationMath.anharmonicBasis(N: 40)
        for (lam, pertTruth, numTruth) in [(0.05, 0.530937, 0.532643),
                                           (0.10, 0.548750, 0.559146)] {
            let e2 = PerturbationVariationMath.secondOrder(state: 0, e0: e0, x4: x4, lam: lam, N: 40)
            let pert = 0.5 + lam * x4[0] + e2
            let num = PerturbationVariationMath.groundEnergy(lam: lam, e0: e0, x4: x4, n: 40)
            #expect(abs(pert - pertTruth) < 5e-7, "λ=\(lam) 微扰 \(pert) vs \(pertTruth)")
            // numTruth 为 Python 打印的 6 位小数锚点，舍入包络 5e-7（绝对）→ rel ~1e-6
            #expect(abs(num - numTruth) / numTruth < 5e-6, "λ=\(lam) 精确 \(num) vs \(numTruth)")
            // 截断微扰在收敛域内逼近精确值（λ 小处相对误差 ~1e-3）
            let rel = abs(num - pert) / num
            #expect(rel < 0.05, "λ=\(lam) 相对误差 \(rel)")
        }
    }

    @Test("物理律：基底收敛——N≥10 的 E₀ 相互一致到 1e-5（截断构造无严格单调，见注释）")
    func basisConvergence() {
        var energies: [Double] = []
        for n in [10, 20, 30, 40] {
            let (e0, x4) = PerturbationVariationMath.anharmonicBasis(N: n)
            energies.append(PerturbationVariationMath.groundEnergy(lam: 0.1, e0: e0, x4: x4, n: n))
        }
        // 注意：不做单调断言——脚本（与移植实现）的 X⁴ 由截断 x 四次连乘（x@x@x@x），
        // 丢掉经由 |N⟩、|N+1⟩ 态的耦合路径，并非 X⁴ 在 span{0..N−1} 上的真投影，
        // 变分上界原理对其不适用：实测 N=10 的 E₀ 反而比收敛值低 2.8e-6（边界路径缺失，
        // 边界附近的 X⁴ 耦合被低估）。有效的物理律是收敛性：N≥10 已一致到 1e-5。
        for e in energies {
            #expect(abs(e - energies[3]) < 1e-5, "E₀=\(e) vs N=40 \(energies[3])")
        }
    }

    // MARK: (A2) 方势阱线性斜坡

    @Test("物理律：一阶微扰高估位移 ~2×（实际位移/一阶 ≈ 0.5，脚本锚点 4.959797）")
    func rampOverestimate() {
        // 精确值锚点（Python 复算）：4.959797
        var h = [Double](repeating: 0, count: 40 * 40)
        for i in 0..<40 {
            h[i * 40 + i] = Double(i + 1) * Double(i + 1) * .pi * .pi / 2.0
                + 0.1 * PerturbationVariationMath.rampMatrixElement(i + 1, i + 1)
            for j in 0..<40 where j != i {
                h[i * 40 + j] = 0.1 * PerturbationVariationMath.rampMatrixElement(i + 1, j + 1)
            }
        }
        let num = AlgebraCore.eigh(h, n: 40).w[0]
        #expect(abs(num - 4.959797) < 1e-5, "精确 E₁ = \(num)")
        let ratio = (num - .pi * .pi / 2.0) / 0.05
        #expect(abs(ratio - 0.5) < 0.01, "位移比 = \(ratio)")
        // 矩阵元解析抽查：x₁₁ = 1/4，x₁₂ = −8/(9π²)，x₁₃ = 0（(1+3) 偶）
        #expect(PerturbationVariationMath.rampMatrixElement(1, 1) == 0.25)
        #expect(abs(PerturbationVariationMath.rampMatrixElement(1, 2) + 8.0 / (9.0 * .pi * .pi)) < 1e-15)
        #expect(PerturbationVariationMath.rampMatrixElement(1, 3) == 0)
    }

    // MARK: (B) 氦变分

    @Test("物理律：上界定理——扫描极小 > 精确非相对论值（−2.903724377 Ha）")
    func variationalUpperBound() {
        let eExact = PerturbationVariationMath.heliumExactHa
        // 网格内全局扫描（400 点，脚本口径）
        var minE = Double.infinity
        for i in 0..<PerturbationVariationMath.variationalScanPoints {
            let a = 1.0 + 1.2 * Double(i) / Double(PerturbationVariationMath.variationalScanPoints - 1)
            minE = min(minE, PerturbationVariationMath.heliumEnergy(a))
        }
        #expect(minE > eExact, "变分极小 \(minE) 应严格大于 \(eExact)")
        // 解析极小：α* = 27/16，E* = (27/16)² − 4·(27/16) + (5/8)(27/16) = −2.84765625（float64 精确）
        let eOpt = PerturbationVariationMath.heliumEnergy(PerturbationVariationMath.heliumAlphaOpt(z: 2.0))
        #expect(abs(eOpt + 2.84765625) < 1e-12)
        for alpha in stride(from: 1.0, through: 2.2, by: 0.05) {
            #expect(PerturbationVariationMath.heliumEnergy(alpha) >= eOpt - 1e-12)
        }
        // 相对偏差 |E*−E_exact|/|E_exact| = 1.93125%（脚本口径取绝对值）
        let rel = abs(eOpt - eExact) / abs(eExact)
        #expect(abs(rel - 0.0193125) < 1e-5, "相对偏差 = \(rel)")
    }

    @Test("W13i 物理律：类氦 Z 扫描——四个离子的变分极小均压住精确值（上界定理逐 Z 成立）")
    func heliumLikeUpperBoundPerZ() {
        // 精确值（Pekeris 型，无限核质量）：H⁻/He/Li⁺/Be²⁺
        let exact: [Double: Double] = [1: -0.5277510165, 2: -2.9037243770,
                                       3: -7.2799134127, 4: -13.6555662400]
        for (z, ex) in exact {
            let aOpt = PerturbationVariationMath.heliumAlphaOpt(z: z)
            let eOpt = PerturbationVariationMath.heliumEnergy(aOpt, z: z)
            #expect(abs(aOpt - (z - 5.0 / 16.0)) < 1e-12)
            #expect(abs(eOpt - (-(z - 5.0 / 16.0) * (z - 5.0 / 16.0))) < 1e-12,
                    "E*(Z=\(z)) = \(eOpt)")
            #expect(eOpt > ex, "Z=\(z)：变分 \(eOpt) 应 > 精确 \(ex)")
            #expect(abs(PerturbationVariationMath.heliumLikeExactHa[z]! - ex) < 1e-9,
                    "模块内建精确表 Z=\(z) 与测试锚点不一致")
        }
    }

    // MARK: fixtures 对拍

    @Test("fixtures：氦 E(α) 扫描 400 点逐点（rel < 1e-8，闭式曲线 + 9 位存储量化）")
    func heliumFixture() throws {
        let fx = try ModuleFixture.load("微扰论与变分法_方势阱与氦变分__v2022")
        let (x, y) = try fx.line(1, 0, index: 0)   // "$E(\alpha)=..."
        #expect(x.count == 400)
        // 闭式曲线本可 1e-12，但 fixture 只存 9 位有效数字（量化 ~3e-9），按标准口径 1e-8
        for i in x.indices {
            let e = PerturbationVariationMath.heliumEnergy(x[i])
            expectPointwiseClose([e], [y[i]], tolerance: 1e-8, "E(α=\(x[i]))")
        }
        // 变分极小标记点（1 点）
        let (mx, my) = try fx.line(1, 0, index: 1)
        #expect(mx.count == 1)
        #expect(abs(mx[0] - 27.0 / 16.0) < 1e-12)
        #expect(abs(my[0] - PerturbationVariationMath.heliumEnergy(27.0 / 16.0)) < 1e-6)
    }

    // MARK: compute 结构与预算

    @Test("compute 输出结构：3 图（误差扫描+λ线；能量级数+λ线，W13g；氦扫描+3 参考线）+ 摘要 5 条")
    func computeStructure() async throws {
        let module = PerturbationVariationModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 3)
        guard case .lineSeries(let c0) = result.charts[0],
              case .lineSeries(let cE) = result.charts[1],
              case .lineSeries(let c1) = result.charts[2] else {
            Issue.record("应为 3 lineSeries"); return
        }
        #expect(c0.series.count == 2, "误差扫描 + 当前 λ 点")
        #expect(c0.series[0].points.count == 30)
        #expect(c0.series[1].points.count == 1)
        #expect(abs(c0.series[1].points[0].x - 0.1) < 1e-12, "默认 λ = 0.1")
        #expect(c0.referenceLines.count == 1 && c0.referenceLines[0].axis == .x)
        // W13g #19：能量级数图（随 λ 联动）
        #expect(cE.series.count == 3)
        #expect(cE.series.allSatisfy { $0.points.count == 30 })
        #expect(cE.series[0].name.contains("一阶") && cE.series[2].name.contains("精确"))
        #expect(cE.referenceLines.count == 1 && cE.referenceLines[0].axis == .x)
        #expect(c1.series[0].points.count == 400)
        #expect(c1.referenceLines.count == 1)
        #expect(abs(c1.referenceLines[0].value - (-2.9037243770)) < 1e-8, "He 精确线")
        #expect(c1.pointMarkers.count == 1)
        #expect(abs(c1.pointMarkers[0].x - 27.0 / 16.0) < 1e-12, "α* = 1.6875")
        #expect(abs(c1.pointMarkers[0].y - (-2.84765625)) < 1e-6, "E* = −2.847656")
        #expect(result.summary.count == 5)
        #expect(result.summary[3].value.contains("满足"), "上界定理成立")
        #expect(result.theory?.formulas.count == 5)
    }

    @Test("W13i #19：Z 联动——Be²⁺（Z=4）扫描窗移到 α*±0.7、极小点随迁、上界仍成立")
    func computeHeliumZLinkage() async throws {
        let module = PerturbationVariationModule()
        var values = ParamValues.defaults(for: module.params)
        values.discretes["heliumZ"] = "4"
        let result = try await module.compute(values, constants: try constants())
        guard case .lineSeries(let c1) = result.charts[2] else {
            Issue.record("氦图应为 lineSeries"); return
        }
        // 扫描窗 [α*−0.7, α*+0.7] = [2.9875, 4.3875]，400 点
        let scan = c1.series[0].points
        #expect(scan.count == 400)
        #expect(abs(scan[0].x - 2.9875) < 1e-12 && abs(scan[399].x - 4.3875) < 1e-12)
        // 极小点：α* = 3.6875，E* = −13.597656
        #expect(abs(c1.pointMarkers[0].x - 3.6875) < 1e-12)
        #expect(abs(c1.pointMarkers[0].y - (-13.59765625)) < 1e-6)
        // 精确线换成 Be²⁺ 值
        #expect(abs(c1.referenceLines[0].value - (-13.6555662400)) < 1e-8)
        // 扫描极小 > 精确（上界）
        let minScan = scan.map { $0.y }.min()!
        #expect(minScan > -13.6555662400)
        #expect(result.summary[2].title.contains("Z=4"))
        #expect(result.summary[3].value.contains("满足"))
    }

    @Test("W13g 物理律：能量级数排序——E⁰¹² < E_exact < E¹ 全扫描域成立（E²<0，一阶高估）")
    func energySeriesOrdering() async throws {
        let module = PerturbationVariationModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        guard case .lineSeries(let cE) = result.charts[1] else {
            Issue.record("图 2 应为 lineSeries"); return
        }
        let e1 = cE.series[0].points, p2 = cE.series[1].points, ex = cE.series[2].points
        for i in e1.indices {
            // E² = λ²Σ|X⁴₀ₘ|²/(E⁰₀−E⁰ₘ) < 0 ⇒ 二阶截断低估；一阶 0.5+λ⟨X⁴⟩₀ 高估
            // （基态为压低 ⟨X⁴⟩ 而展宽）；精确值居中。锚点 λ=0.1：0.5488 < 0.5591 < 0.5750。
            #expect(p2[i].y < ex[i].y && ex[i].y < e1[i].y,
                    "λ=\(e1[i].x)：E⁰¹² \(p2[i].y) / exact \(ex[i].y) / E¹ \(e1[i].y)")
        }
    }

    @Test("秒级档预算：compute < 2000 ms（31 次 40² eigh + 扫描，强制口径）")
    func computeBudget() async throws {
        let module = PerturbationVariationModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 2000)
    }
}
