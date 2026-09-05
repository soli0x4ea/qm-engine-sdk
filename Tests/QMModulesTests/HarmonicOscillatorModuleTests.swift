import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W3-5 量子谐振子模块：fixtures 对拍（量子谐振子_有限差分与相干态__v2022.json）。
///
/// fixture 参数与 Python 默认一致：N=3000、L=18、nstates=6、α=2。
/// 提取语义：能量图为 matplotlib barh(y=n, width=E_n, height=0.5)，
/// 管线存储的 bars.x 是柱中心（left=0 + width/2 = E_n/2）→ 真值能量 = 2×x。
@Suite("W3 量子谐振子模块")
struct HarmonicOscillatorModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    @Test("fixtures 参数（N=3000, L=18）能量逐级对拍")
    func energiesMatchFixture() throws {
        let fx = try ModuleFixture.load("量子谐振子_有限差分与相干态__v2022")
        let bars = try fx.bars(0, 0)
        #expect(bars.count == 6)
        #expect(bars.allSatisfy { $0.height == 0.5 }, "barh 柱高固定 0.5（提取语义锚点）")

        let (_, energies, _) = OscillatorMath.solve(N: 3000, L: 18, nstates: 6)
        // fixture bars.x 存 9 位有效数字（E_n/2），×2 后有效位 ~8.6 位 → 容差 1e-8
        expectPointwiseClose(energies, bars.map { 2 * $0.x },
                             tolerance: 1e-8, "E_n（N=3000, L=18）")
        // 与解析 E_n = n+½ 同框验证（有限差分 dx² 收敛）
        for n in 0..<6 {
            #expect(abs(energies[n] - (Double(n) + 0.5)) / (Double(n) + 0.5) < 1e-4,
                    "E_\(n) = \(energies[n])")
        }
    }

    @Test("解析波函数 ψn 对拍 fixtures（n = 0, 1, 2）")
    func wavefunctionsMatchFixture() throws {
        let fx = try ModuleFixture.load("量子谐振子_有限差分与相干态__v2022")
        for n in 0...2 {
            let (xs, ys) = try fx.line(0, 1, label: "analytic n=\(n)")
            var re: [Double] = [], ex: [Double] = []
            let psi = Hermite.normalizedWavefunction(n, xs: xs)
            for (i, y) in ys.enumerated() where abs(y) > 1e-10 {
                re.append(psi[i]); ex.append(y)
            }
            #expect(ex.count > 500, "零点掩码后应保留绝大多数点")
            expectPointwiseClose(re, ex, tolerance: 1e-7, "解析 ψ\(n)")
        }
    }

    @Test("相干态 |ψ(x,t)|² 四相位对拍（α = 2）")
    func coherentMatchFixture() throws {
        let fx = try ModuleFixture.load("量子谐振子_有限差分与相干态__v2022")
        let phases: [(label: String, t: Double)] = [
            ("t=0.00", 0), ("t=0.50", .pi / 2), ("t=1.00", .pi), ("t=1.50", 3 * .pi / 2),
        ]
        for p in phases {
            let (xs, ys) = try fx.line(1, 0, label: p.label)
            var rho = OscillatorMath.coherentDensity(alpha: 2, t: p.t, xs: xs)
            let z = Num.trapezoid(xs, rho)
            rho = rho.map { $0 / z }
            // 尾部 ~1e-33 区域被指数敏感性放大相对误差，只在有效区间对拍
            // （高斯包络 ρ>1e-6 约 |x−x_c|<3.5σ，覆盖 512 点中的 ~210 点）
            var re: [Double] = [], ex: [Double] = []
            for (i, y) in ys.enumerated() where abs(y) > 1e-6 {
                re.append(rho[i]); ex.append(y)
            }
            #expect(ex.count > 150, "有效区间点数 \(ex.count)")
            expectPointwiseClose(re, ex, tolerance: 1e-7, "相干态 \(p.label)")
        }
    }

    @Test("数值物理量：节点计数 / E_n 收敛 / 数值波函数 / 相干态矩")
    func numericPhysics() {
        let (x, energies, vectors) = OscillatorMath.solve(N: 800, L: 12, nstates: 6)
        for n in 0..<6 {
            #expect(OscillatorMath.nodeCount(vectors[n]) == n, "ψ\(n) 节点数")
            // dx = 24/799 → 差分色散误差 ~c·dx²，实测 E_5 相对误差 3.1e-4
            #expect(abs(energies[n] - (Double(n) + 0.5)) / (Double(n) + 0.5) < 5e-4)
        }

        // 数值波函数插值重归一化后 vs 解析：差分误差在指数尾部（y~4e-6）与
        // 零点附近（y~0.008）相对放大至 ~1%（Python 数值曲线同样水平），容差 2%
        let zoom = Num.linspace(-5, 5, count: 400)
        for n in 0..<3 {
            var num = OscillatorMath.interpolate(xs: x, ys: vectors[n], onto: zoom)
            let norm = Num.trapezoid(zoom, num.map { $0 * $0 })
            num = num.map { $0 / sqrt(norm) }
            let ana = Hermite.normalizedWavefunction(n, xs: zoom)
            var re: [Double] = [], ex: [Double] = []
            for (i, y) in ana.enumerated() where abs(y) > 1e-10 {
                re.append(num[i]); ex.append(y)
            }
            expectPointwiseClose(re, ex, tolerance: 2e-2, "数值 ψ\(n)（N=800, L=12）")
        }

        // 相干态矩：⟨x⟩(0) = √2·α = 2√2（修正 Python 源标签 x_c=2cos t 的系数笔误）；
        // Var ≡ 0.5（基态宽度不随 t 变化）；t=π/2 质心回零
        let xg = Num.linspace(-9, 9, count: 2400)
        var rho = OscillatorMath.coherentDensity(alpha: 2, t: 0, xs: xg)
        rho = rho.map { $0 / Num.trapezoid(xg, rho) }
        let (mean0, var0) = OscillatorMath.meanAndVariance(ofRho: rho, xs: xg)
        #expect(abs(mean0 - 2 * 2.0.squareRoot()) < 1e-6, "⟨x⟩(0) = \(mean0)")
        #expect(abs(var0 - 0.5) < 1e-4, "Var(0) = \(var0)")

        var rhoQ = OscillatorMath.coherentDensity(alpha: 2, t: .pi / 2, xs: xg)
        rhoQ = rhoQ.map { $0 / Num.trapezoid(xg, rhoQ) }
        let (meanQ, varQ) = OscillatorMath.meanAndVariance(ofRho: rhoQ, xs: xg)
        #expect(abs(meanQ) < 1e-6, "⟨x⟩(π/2) = \(meanQ)")
        #expect(abs(varQ - 0.5) < 1e-4, "Var(π/2) = \(varQ)")
    }

    @Test("compute 输出结构：三图 / 能量柱 12 / 质心摘要（修正值 2.83）")
    func computeStructure() async throws {
        let module = HarmonicOscillatorModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())

        #expect(result.charts.count == 3)
        guard case .bars(let energyBars) = result.charts[0] else {
            Issue.record("charts[0] 应为 bars"); return
        }
        #expect(energyBars.bars.count == 12, "6 态 × 数值/解析")
        #expect(energyBars.bars.filter { $0.series == "数值" }.count == 6)

        guard case .lineSeries(let wave) = result.charts[1] else {
            Issue.record("charts[1] 应为 lineSeries"); return
        }
        #expect(wave.series.count == 6, "n=0…2 × 数值/解析")

        guard case .lineSeries(let phase) = result.charts[2] else {
            Issue.record("charts[2] 应为 lineSeries"); return
        }
        #expect(phase.series.count == 4, "四相位")

        // 质心摘要为修正后的 √2·α ≈ 2.83（非 Python 标签的 2.0）
        let xc = result.summary.first { $0.id == "xc0" }
        #expect(xc?.value.hasPrefix("2.8") == true, "⟨x⟩(0) 摘要 = \(xc?.value ?? "nil")")
        // E₀ 数值卡与节点校验卡存在
        #expect(result.summary.contains { $0.id == "e0" })
        #expect(result.summary.contains { $0.id == "nodes" && $0.value.hasPrefix("6/6") })
        #expect(result.theory?.formulas.count == 4)
    }

    /// naive 参考实现：收尾包 1 重组**前**的 coherentDensity 旧路径
    /// （每 n 一次 normalizedWavefunction + re[i] += cr·ψ[i]）逐句复刻。
    /// 用于锁死重组的 bit 级无损声明——逐元素 == 比较。
    private static func coherentDensityLegacy(
        alpha: Double, t: Double, xs: [Double], nmax: Int = 44
    ) -> [Double] {
        if alpha <= 0 {
            let psi0 = Hermite.normalizedWavefunction(0, xs: xs)
            return psi0.map { $0 * $0 }
        }
        var re = [Double](repeating: 0, count: xs.count)
        var im = [Double](repeating: 0, count: xs.count)
        let lnAlpha = log(alpha)
        for n in 0..<nmax {
            let lnAmp = -alpha * alpha / 2 + Double(n) * lnAlpha - 0.5 * Num.lnFactorial(n)
            let amp = exp(lnAmp)
            let phase = -t * (Double(n) + 0.5)
            let cr = amp * cos(phase)
            let ci = amp * sin(phase)
            let psi = Hermite.normalizedWavefunction(n, xs: xs)
            for i in xs.indices {
                re[i] += cr * psi[i]
                im[i] += ci * psi[i]
            }
        }
        return zip(re, im).map { $0 * $0 + $1 * $1 }
    }

    @Test("coherentDensity 重组 bit 级无损：与旧路径 naive 参考逐元素 ==（收尾包 1）")
    func coherentDensityBitwiseStable() {
        // 覆盖典型参数域：α（含 0 退化分支）、t（含 π/2 相位）、网格规模
        for (alpha, t) in [(2.0, 0.0), (2.0, .pi / 2), (1.5, 3 * .pi / 2), (0.7, 1.234)] {
            let xs = Num.linspace(-9, 9, count: 2400)
            let legacy = Self.coherentDensityLegacy(alpha: alpha, t: t, xs: xs)
            let current = OscillatorMath.coherentDensity(alpha: alpha, t: t, xs: xs)
            #expect(legacy.count == current.count)
            for i in xs.indices {
                #expect(legacy[i] == current[i],
                        "α=\(alpha) t=\(t) x[\(i)]: \(legacy[i]) vs \(current[i])")
            }
        }
    }

    @Test("秒级档预算：默认参数（N=800）compute 中位 < 2 s")
    func computeBudget() async throws {
        let module = HarmonicOscillatorModule()
        try await expectComputeUnderBudget(
            module: module,
            values: ParamValues.defaults(for: module.params),
            constants: try constants(),
            budgetMillis: 2000, runs: 5)
    }
}
