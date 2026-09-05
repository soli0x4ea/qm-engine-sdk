import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W4-5 高斯波包演化：fixtures 曲线对拍（512 点 ×8 曲线）+ 脚本断言复刻。
@Suite("W4 高斯波包演化")
struct WavePacketModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }
    private static let sigma0 = 0.2e-6

    @Test("check1：v_g = v、v_p = v/2（rel < 1e-12）")
    func velocities() throws {
        let cs = try constants()
        let hbar = try cs.value("hbar"), me = try cs.value("m_e"), c = try cs.value("c")
        let v = 0.01 * c
        let k0 = me * v / hbar
        let vg = WavePacketMath.groupVelocity(k0: k0, hbar: hbar, m: me)
        let vp = WavePacketMath.phaseVelocity(k0: k0, hbar: hbar, m: me)
        #expect(abs(vg - v) / v < 1e-12, "v_g=\(vg)")
        #expect(abs(vp - v / 2) / (v / 2) < 1e-12, "v_p=\(vp)")
        // 德布罗意波长量级（0.01c → ~242.6 pm）
        #expect(abs(2 * .pi / k0 * 1e12 - 242.63) < 0.5)
    }

    @Test("check2：实验室系峰值位置 = v_g·t（abs < 1e-8）且宽度不变（rel < 1e-3）")
    func labFrameTranslation() throws {
        let cs = try constants()
        let hbar = try cs.value("hbar"), me = try cs.value("m_e"), c = try cs.value("c")
        let vg = 0.01 * c
        let xlab = Num.linspace(-1e-6, 9e-6, count: 2000)
        for tps in WavePacketModule.labTimesPs {
            let t = tps * 1e-12
            let xc = vg * t
            let st = WavePacketMath.sigmaT(t, sigma0: Self.sigma0, hbar: hbar, m: me)
            let p = xlab.map { WavePacketMath.psi2(x: $0, xc: xc, sigma: st) }
            let iPeak = p.indices.max(by: { p[$0] < p[$1] })!
            let absPos = abs(xlab[iPeak] - xc)
            #expect(absPos < 1e-8, "t=\(tps)ps 峰值偏差 \(absPos)")
            let sigNum = WavePacketMath.numericWidth(xs: xlab, psi2: p, xpeak: xlab[iPeak])
            #expect(abs(sigNum - st) / st < 1e-3, "t=\(tps)ps 宽度 rel")
        }
    }

    @Test("check3：共动系数值宽度 = σ(t)（4 个时刻 rel < 1e-3）")
    func comovingSpreading() throws {
        let cs = try constants()
        let hbar = try cs.value("hbar"), me = try cs.value("m_e")
        let xp = Num.linspace(-1e-6, 1.5e-6, count: 2000)
        for tns in WavePacketModule.comTimesNs {
            let t = tns * 1e-9
            let st = WavePacketMath.sigmaT(t, sigma0: Self.sigma0, hbar: hbar, m: me)
            let p = xp.map { WavePacketMath.psi2(x: $0, xc: 0, sigma: st) }
            let sigNum = WavePacketMath.numericWidth(xs: xp, psi2: p, xpeak: 0)
            #expect(abs(sigNum - st) / st < 1e-3, "t=\(tns)ns rel=\(abs(sigNum - st) / st)")
        }
    }

    @Test("fixture 曲线对拍：实验室系 4 条 + 共动系 4 条（512 点）")
    func curvesMatchFixtures() throws {
        let fx = try ModuleFixture.load("物质波包_高斯波包演化__v2022")
        let cs = try constants()
        let hbar = try cs.value("hbar"), me = try cs.value("m_e"), c = try cs.value("c")
        let vg = 0.01 * c

        // fig0ax0：实验室系 4 曲线（x 单位 μm，y 归一化）
        for (i, tps) in WavePacketModule.labTimesPs.enumerated() {
            let (fx_, fy) = try fx.line(0, 0, index: i)
            let t = tps * 1e-12
            let st = WavePacketMath.sigmaT(t, sigma0: Self.sigma0, hbar: hbar, m: me)
            // 归一化口径与脚本一致：除以该曲线自身峰值（2000 点网格）
            let xc = vg * t
            let grid = Num.linspace(-1e-6, 9e-6, count: 2000)
            let raw = grid.map { WavePacketMath.psi2(x: $0, xc: xc, sigma: st) }
            let mx = raw.max() ?? 1
            expectPointwiseClose(
                fx_.map { WavePacketMath.psi2(x: $0 * 1e-6, xc: xc, sigma: st) / mx },
                fy, tolerance: 1e-6, "实验室系 t=\(tps)ps")
        }
        // fig0ax1：共动系 4 曲线
        for (i, tns) in WavePacketModule.comTimesNs.enumerated() {
            let (fx_, fy) = try fx.line(0, 1, index: i)
            let t = tns * 1e-9
            let st = WavePacketMath.sigmaT(t, sigma0: Self.sigma0, hbar: hbar, m: me)
            let grid = Num.linspace(-1e-6, 1.5e-6, count: 2000)
            let raw = grid.map { WavePacketMath.psi2(x: $0, xc: 0, sigma: st) }
            let mx = raw.max() ?? 1
            expectPointwiseClose(
                fx_.map { WavePacketMath.psi2(x: $0 * 1e-6, xc: 0, sigma: st) / mx },
                fy, tolerance: 1e-6, "共动系 t=\(tns)ns")
        }
    }

    @Test("compute 输出结构：双图组 4+4 系列 + 摘要 + 理论卡")
    func computeStructure() async throws {
        let module = WavePacketModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        guard case .lineSeries(let lab) = result.charts[0],
              case .lineSeries(let com) = result.charts[1] else { return }
        #expect(lab.series.count == 4)
        #expect(com.series.count == 4)
        #expect(lab.series[0].points.count == 500)   // 2000 / 4
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("实时档预算：compute 中位 < 16 ms")
    func computeBudget() async throws {
        let module = WavePacketModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }
}
