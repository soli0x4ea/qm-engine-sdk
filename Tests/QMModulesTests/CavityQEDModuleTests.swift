import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W6 笔记 36《腔量子电动力学与量子比特实现》：JC 真空 Rabi 振荡 / 受驱量子比特退相干。
/// fixtures：腔量子电动力学_JaynesCummings__v2022.json（P_e 双平台 512 点 + 归一化吸收谱）、
/// 腔量子电动力学_量子比特Rabi__v2022.json（三平台 512 点）。
@Suite("W6 腔量子电动力学")
struct CavityQEDModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    // 脚本固定口径（rad/s）
    private static let gCQED = 2.0 * .pi * 47e3
    private static let gCircuit = 2.0 * .pi * 5e6
    private static let omegaRDrive = 2.0 * .pi * 1e6

    // MARK: - 物理律

    @Test("物理律：真空 Rabi 振荡 P_e(0)=1、半周期完全转移、周期 T_R = π/g")
    func vacuumRabiLaw() {
        let g = Self.gCQED
        #expect(abs(CavityQEDMath.vacuumRabiPe(0, g: g) - 1) < 1e-15, "初始全激发")
        #expect(abs(CavityQEDMath.vacuumRabiPe(.pi / (2 * g), g: g)) < 1e-15,
                "t = π/2g 处能量完全转移到腔（P_e = 0）")
        #expect(abs(CavityQEDMath.vacuumRabiPe(.pi / g, g: g) - 1) < 1e-15,
                "t = π/g 处回到原子（周期 T_R = π/g）")
        #expect(abs(CavityQEDMath.oscillationPeriod(g: g) - .pi / g) < 1e-18)
        // 真空 Rabi 频率 Ω_R = 2g；P_e = (1 + cos Ω_R t)/2 恒等
        for t in stride(from: 0.0, to: 1e-4, by: 1.3e-6) {
            #expect(abs(CavityQEDMath.vacuumRabiPe(t, g: g)
                        - (1 + cos(2 * g * t)) / 2) < 1e-15)
        }
    }

    @Test("物理律：JC Rabi 频率 √(n+1) 标度——Ω_1/Ω_0 = √2、Ω_4/Ω_0 = √5")
    func rabiScalingLaw() {
        let g = Self.gCircuit
        #expect(abs(CavityQEDMath.rabiFrequency(g: g, n: 1)
                    / CavityQEDMath.rabiFrequency(g: g, n: 0) - 2.0.squareRoot()) < 1e-15)
        #expect(abs(CavityQEDMath.rabiFrequency(g: g, n: 4)
                    / CavityQEDMath.rabiFrequency(g: g, n: 0) - 5.0.squareRoot()) < 1e-15)
        #expect(abs(CavityQEDMath.rabiFrequency(g: g, n: 0) - 2 * g) < 1e-9, "真空 Ω_R = 2g")
    }

    @Test("物理律：缀饰态吸收谱双峰在 δ = ±g，峰距 = 2g（强耦合正常模分裂）")
    func dressedSplittingLaw() {
        let g = Self.gCircuit
        let kappa = 2.0 * .pi * 0.5e6
        // 网格上找峰位：双峰对称于 δ = 0
        let step = 1.6 * g / 400
        let delta = Array(stride(from: -1.6 * g, through: 1.6 * g, by: step))
        let a = delta.map { CavityQEDMath.absorption($0, g: g, kappa: kappa) }
        let peakValue = a.max() ?? 0
        let peakDelta = delta[a.firstIndex(of: peakValue) ?? 0]
        #expect(abs(abs(peakDelta) - g) < 2 * step,
                "峰位 |δ*| = \(abs(peakDelta)) ≈ g = \(g)（网格分辨内）")
        // 对称性 A(−δ) = A(+δ)
        for d in stride(from: 0.1 * g, through: 1.5 * g, by: 0.13 * g) {
            #expect(abs(CavityQEDMath.absorption(d, g: g, kappa: kappa)
                        - CavityQEDMath.absorption(-d, g: g, kappa: kappa)) < 1e-6)
        }
        // 共振点为深谷：A(0) ≪ A(g)（分辨的双峰）
        let a0 = CavityQEDMath.absorption(0, g: g, kappa: kappa)
        let ag = CavityQEDMath.absorption(g, g: g, kappa: kappa)
        #expect(ag / a0 > 10, "g/κ = 10：强耦合双峰可分辨（比值 \(ag / a0)）")
        // 洛伦兹主项自检：A(g) ≈ g²/(κ/2)²（另一峰尾贡献 ~1/4，相对 < 1e-3）
        #expect(abs(ag - g * g / (kappa / 2 * (kappa / 2))) < 1e-3 * ag)
    }

    @Test("物理律：退相干包络上界 P_e ≤ e^{−t/T₂}；首峰恰在 t = π/Ω_R")
    func decoherenceEnvelopeLaw() {
        let omegaR = Self.omegaRDrive
        let t2 = 50e-6   // transmon
        for t in stride(from: 0.0, to: 200e-6, by: 7.7e-6) {
            let pe = CavityQEDMath.drivenQubitPe(t, omegaR: omegaR, t2: t2)
            #expect(pe <= exp(-t / t2) + 1e-15, "t = \(t)：sin² ≤ 1 给出包络")
            #expect(pe >= -1e-15, "概率非负")
        }
        #expect(abs(CavityQEDMath.drivenQubitPe(0, omegaR: omegaR, t2: t2)) < 1e-15,
                "t = 0 从基态出发")
        // 首峰：sin²(Ω t/2) = 1 ⟺ t = π/Ω；P_e = e^{−π/(Ω T₂)}
        let tPeak = Double.pi / omegaR
        let expected = exp(-tPeak / t2)
        #expect(abs(CavityQEDMath.drivenQubitPe(tPeak, omegaR: omegaR, t2: t2) - expected) < 1e-15)
        // 相干振荡次数 N_osc = Ω_R T₂/π：离子阱 T₂=1 s → 2×10⁶ 次
        #expect(abs(CavityQEDMath.coherentOscillationCount(omegaR: omegaR, t2: 1.0)
                    - 2.0e6) < 1e-9)
        // 平台相干时间排序：离子阱 > NV > transmon（计算预算的量级差异）
        let ion = QubitRabiModule.platforms[0].t2
        let sc = QubitRabiModule.platforms[1].t2
        let nv = QubitRabiModule.platforms[2].t2
        #expect(ion > nv && nv > sc)
    }

    // MARK: - fixtures 对拍

    @Test("fixtures：JC 真空 Rabi P_e 双平台 512 点（混合容差：|Δ| ≤ max(1e-9, 2e-5·|y|)）")
    func vacuumRabiMatchesFixtures() throws {
        let fx = try ModuleFixture.load("腔量子电动力学_JaynesCummings__v2022")
        // fixtures x 为 μs 且 9 位有效数字存储——t 的舍入经 g 放大后产生 ~1e-5 相对差
        let (tCav, yCav) = try fx.line(0, 0, label: "cavity QED")
        let (tCir, yCir) = try fx.line(0, 0, label: "circuit QED")
        #expect(tCav.count == 512 && tCir.count == 512)
        for i in 0..<512 {
            let p = CavityQEDMath.vacuumRabiPe(tCav[i] * 1e-6, g: Self.gCQED)
            #expect(abs(p - yCav[i]) <= max(1e-9, 2e-5 * abs(yCav[i])),
                    "cavity 第 \(i) 点：\(p) vs \(yCav[i])")
            let q = CavityQEDMath.vacuumRabiPe(tCir[i] * 1e-6, g: Self.gCircuit)
            #expect(abs(q - yCir[i]) <= max(1e-9, 2e-5 * abs(yCir[i])),
                    "circuit 第 \(i) 点：\(q) vs \(yCir[i])")
        }
    }

    @Test("fixtures：缀饰态吸收谱归一化曲线 512 点 ×2（< 1e-7，按 fixture 网格自身峰值归一）")
    func absorptionMatchesFixtures() throws {
        let fx = try ModuleFixture.load("腔量子电动力学_JaynesCummings__v2022")
        // 电路 QED：x 轴 MHz，κ = 2π·0.5 MHz
        let (xCir, yCir) = try fx.line(0, 1, label: "circuit QED")
        #expect(xCir.count == 512)
        let aCir = xCir.map { CavityQEDMath.absorption($0 * 2 * .pi * 1e6, g: Self.gCircuit,
                                                        kappa: 2 * .pi * 0.5e6) }
        let maxCir = aCir.max() ?? 1
        for i in 0..<512 {
            #expect(abs(aCir[i] / maxCir - yCir[i]) < 1e-7,
                    "circuit 谱 第 \(i) 点：\(aCir[i] / maxCir) vs \(yCir[i])")
        }
        // 腔 QED：x 轴 kHz，κ = 2π·5 kHz
        let (xCav, yCav) = try fx.line(0, 1, label: "cavity QED")
        #expect(xCav.count == 512)
        let aCav = xCav.map { CavityQEDMath.absorption($0 * 2 * .pi * 1e3, g: Self.gCQED,
                                                        kappa: 2 * .pi * 5e3) }
        let maxCav = aCav.max() ?? 1
        for i in 0..<512 {
            #expect(abs(aCav[i] / maxCav - yCav[i]) < 1e-7,
                    "cavity 谱 第 \(i) 点：\(aCav[i] / maxCav) vs \(yCav[i])")
        }
    }

    @Test("fixtures：受驱量子比特三平台 512 点（|Δ| < 2e-5，t 的 9 位存储舍入经 Ω_R 放大）")
    func qubitRabiMatchesFixtures() throws {
        let fx = try ModuleFixture.load("腔量子电动力学_量子比特Rabi__v2022")
        let t2s = [("Ion trap", 1.0), ("Superconducting", 50e-6), ("NV center", 1e-3)]
        for (label, t2) in t2s {
            let (t, y) = try fx.line(0, 0, label: label)
            #expect(t.count == 512)
            for i in 0..<512 {
                let p = CavityQEDMath.drivenQubitPe(t[i] * 1e-6, omegaR: Self.omegaRDrive, t2: t2)
                #expect(abs(p - y[i]) < 2e-5, "\(label) 第 \(i) 点：\(p) vs \(y[i])")
            }
        }
    }

    // MARK: - compute 结构与预算

    @Test("compute 结构：JC 双图（600 点双平台 + 800 点双谱）+ 摘要 + 理论卡；实时档预算")
    func jcComputeStructure() async throws {
        let module = JaynesCummingsModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        for (index, chart) in result.charts.enumerated() {
            guard case .lineSeries(let c) = chart else {
                Issue.record("charts[\(index)] 应为 lineSeries"); return
            }
            #expect(c.series.count == 2)
            #expect(c.series[0].points.count == (index == 0 ? 600 : 800))
        }
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 4)
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }

    @Test("compute 结构：量子比特 Rabi 三平台（各 500 点）+ 摘要 + 理论卡；实时档预算")
    func qubitRabiComputeStructure() async throws {
        let module = QubitRabiModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 1)
        guard case .lineSeries(let c) = result.charts[0] else {
            Issue.record("应为 lineSeries"); return
        }
        #expect(c.series.count == 3)
        for s in c.series { #expect(s.points.count == 500, "2000/4 抽稀") }
        #expect(result.summary.count == 5)
        #expect(result.theory?.formulas.count == 4)
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }
}
