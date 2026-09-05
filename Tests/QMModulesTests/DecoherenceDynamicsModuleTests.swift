import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W5-13 退相干动力学：Caldeira-Leggett τ_D 尺度 + T2=2T1 / 指针基 / 猫态负性
/// （fixtures 4 图 10 条曲线全量 400 点对拍；脚本无 [check]，公式复算口径）。
@Suite("W5 量子测量：退相干动力学")
struct DecoherenceDynamicsModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    @Test("τ_D 尺度：1 nm → 3.957e-18 s（默认 η、T）；a⁻³ 律每量级 ×10⁻³")
    func tauD() throws {
        let cs = try constants()
        let hbar = try cs.value("hbar"), kB = try cs.value("kB")
        let t9 = DecoherenceDynamicsMath.tauD(1e-9, eta: 1.8e-5, T: 300, hbar: hbar, kB: kB)
        #expect(abs(t9 / 3.9568e-18 - 1) < 1e-3, "τ_D(1 nm) = \(t9)")
        let t6 = DecoherenceDynamicsMath.tauD(1e-6, eta: 1.8e-5, T: 300, hbar: hbar, kB: kB)
        let t3 = DecoherenceDynamicsMath.tauD(1e-3, eta: 1.8e-5, T: 300, hbar: hbar, kB: kB)
        #expect(abs(t6 / t9 - 1e-9) < 1e-15, "a⁻³：1 µm 是 1 nm 的 10⁻⁹")
        #expect(abs(t3 / t6 - 1e-9) < 1e-15)
        // 环境温度升高 2 倍 → 退相干加速 2 倍（τ_D ∝ 1/T）
        let t300 = DecoherenceDynamicsMath.tauD(1e-9, eta: 1.8e-5, T: 300, hbar: hbar, kB: kB)
        let t600 = DecoherenceDynamicsMath.tauD(1e-9, eta: 1.8e-5, T: 600, hbar: hbar, kB: kB)
        #expect(abs(t300 / t600 - 2) < 1e-12)
    }

    @Test("两能级：|ρ_eg(t)|² = ρ_ee(t)（T2 = 2T1）")
    func twoLevel() {
        let t = [0.0, 0.5, 1.0, 2.5, 5.0]
        for ti in t {
            let ee = DecoherenceDynamicsMath.rhoEE(ti, gamma: 1)
            let eg = DecoherenceDynamicsMath.rhoEG(ti, gamma: 1)
            #expect(abs(eg * eg - ee) < 1e-14, "t = \(ti)")
        }
        #expect(DecoherenceDynamicsMath.rhoEE(0, gamma: 1) == 1)
        // Γ 翻倍 → 指数律 e^(−2γt) = (e^(−γt))²；相干元 Γt 组合不变
        #expect(abs(DecoherenceDynamicsMath.rhoEE(1, gamma: 2)
                - pow(DecoherenceDynamicsMath.rhoEE(1, gamma: 1), 2)) < 1e-15)
        #expect(abs(DecoherenceDynamicsMath.rhoEG(2, gamma: 1)
                - DecoherenceDynamicsMath.rhoEG(1, gamma: 2)) < 1e-15)
    }

    @Test("指针基：叠加相干 e^(−2x₀²t/τ₀)，x₀ = 1、τ₀ = 0.6；对角混态恒为零")
    func pointerBasis() {
        #expect(abs(DecoherenceDynamicsMath.pointerCoherence(0, x0: 1, tau0: 0.6) - 1) < 1e-15)
        let t = 1.0
        #expect(abs(DecoherenceDynamicsMath.pointerCoherence(t, x0: 1, tau0: 0.6)
                    - exp(-4.0 * t / 1.2)) < 1e-15)
        // 分离越远退相干越快：(2x₀)² 律
        let c1 = DecoherenceDynamicsMath.pointerCoherence(0.5, x0: 1, tau0: 0.6)
        let c2 = DecoherenceDynamicsMath.pointerCoherence(0.5, x0: 2, tau0: 0.6)
        #expect(abs(c2 - pow(c1, 4)) < 1e-14, "(2x₀)² 律")
    }

    @Test("猫态负性：N(0) = 1；Γ_cat = 2κ|α|²，α = 1/2/3 → Γ = 2/8/18")
    func catStates() {
        for (alpha, gamma) in zip([1.0, 2.0, 3.0], [2.0, 8.0, 18.0]) {
            #expect(DecoherenceDynamicsMath.catNegativity(0, kappa: 1, alpha: alpha) == 1)
            let n = DecoherenceDynamicsMath.catNegativity(0.3, kappa: 1, alpha: alpha)
            #expect(abs(n - exp(-gamma * 0.3)) < 1e-14, "α = \(alpha)")
        }
    }

    @Test("fixture 曲线对拍：τ_D(a) loglog 400 点（图 0）")
    func tauCurveFixture() throws {
        let cs = try constants()
        let hbar = try cs.value("hbar"), kB = try cs.value("kB")
        let fx = try ModuleFixture.load("量子测量与退相干_退相干动力学__v2022")
        let (x, y) = try fx.line(0, 0, index: 0)
        #expect(x.count == 400)
        #expect(abs(x.first! - 1e-9) < 1e-20 && abs(x.last! - 1e-3) < 1e-9)
        let tau = x.map { DecoherenceDynamicsMath.tauD($0, eta: 1.8e-5, T: 300,
                                                       hbar: hbar, kB: kB) }
        expectPointwiseClose(tau, y, tolerance: 1e-7, "τ_D(a) = ℏ²/(12πηk_BT·a³)")
        // a³ 律全序列自检：任意点可由首点按 a⁻³ 推出
        for i in stride(from: 1, to: 400, by: 50) {
            let ratio = pow(x[0] / x[i], 3)
            #expect(abs(y[i] / y[0] / ratio - 1) < 1e-6)
        }
    }

    @Test("fixture 曲线对拍：三子图 8 条曲线各 400 点（图 1）")
    func subplotsFixture() throws {
        let fx = try ModuleFixture.load("量子测量与退相干_退相干动力学__v2022")
        // (a) 两能级：ρ_ee = e^(−t)、|ρ_eg| = e^(−t/2)，t ∈ [0, 5]
        let (t, ee) = try fx.line(1, 0, index: 0)
        #expect(t.count == 400)
        expectPointwiseClose(t.map { DecoherenceDynamicsMath.rhoEE($0, gamma: 1) }, ee,
                             tolerance: 1e-7, "ρ_ee(t) = e^(−t/T1)")
        let (_, eg) = try fx.line(1, 0, index: 1)
        expectPointwiseClose(t.map { DecoherenceDynamicsMath.rhoEG($0, gamma: 1) }, eg,
                             tolerance: 1e-7, "|ρ_eg(t)| = e^(−t/2T1)")
        // (b) 指针基：叠加相干衰减 + 混态恒 0，tt ∈ [0, 4]
        let (tt, coh) = try fx.line(1, 1, index: 0)
        #expect(tt.count == 400)
        expectPointwiseClose(tt.map { DecoherenceDynamicsMath.pointerCoherence($0, x0: 1, tau0: 0.6) },
                             coh, tolerance: 1e-7, "指针基叠加相干")
        let (_, mix) = try fx.line(1, 1, index: 1)
        #expect(mix.allSatisfy { $0 == 0 }, "位置对角混态非对角元恒为 0")
        // (c) 猫态负性：α = 1/2/3 三条
        for (idx, alpha) in zip([0, 1, 2], [1.0, 2.0, 3.0]) {
            let (tc, n) = try fx.line(1, 2, index: idx)
            #expect(tc.count == 400)
            expectPointwiseClose(tc.map { DecoherenceDynamicsMath.catNegativity($0, kappa: 1, alpha: alpha) },
                                  n, tolerance: 1e-7, "猫态 |α| = \(alpha)")
        }
    }

    @Test("compute 输出结构：4 图（τ_D loglog + 三子图）+ 摘要 3 条 + 理论卡")
    func computeStructure() async throws {
        let module = DecoherenceDynamicsModule()
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
        #expect(c0.spec.xAxis.scale == .log && c0.spec.yAxis.scale == .log)
        #expect(c0.referenceLines.count == 3, "1 s / 1 ms / 1 ps 时间尺度参考线")
        #expect(c1.series.count == 2 && c1.series[0].points.count == 400)
        #expect(c2.series.count == 2 && c2.series[0].points.count == 400)
        #expect(c3.series.count == 3 && c3.series[0].points.count == 400)
        #expect(result.summary.count == 3)
        #expect(result.summary[1].value.contains("e-27"), "τ_D(1 µm) ~ 1e-27 s：\(result.summary[1].value)")
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("实时档预算：compute 中位 < 16 ms")
    func computeBudget() async throws {
        let module = DecoherenceDynamicsModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }
}
