import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W4-12 算符表象：复刻脚本打印值（Δx/Δp/乘积）+ 闭式高斯对拍 + fixtures 双线对拍。
/// 动量 fixture 的 x 标签存在导出 quirk（中心区域跳格，p=0 峰点丢失），
/// 故对拍口径为「|x| 反推网格下标 + 信号点 y 对拍」，不依赖其 x 标签。
@Suite("W4 算符表象")
struct OperatorRepresentationModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }
    private static let sigmaX = 2.0

    /// 复刻脚本状态：σₓ = 2 的归一化高斯 + FFT 动量表象。
    private static func state() -> (x: [Double], psi2: [Double], p: [Double], phi2: [Double]) {
        let (x, psi, dx) = OperatorRepresentationMath.gaussianState(sigmaX: sigmaX)
        var psi2 = [Double](repeating: 0, count: x.count)
        for i in x.indices { psi2[i] = psi[i] * psi[i] }
        let (p, phi2) = OperatorRepresentationMath.momentumDensity(psi: psi, dx: dx)
        return (x, psi2, p, phi2)
    }

    @Test("脚本打印值：Δx = 1.4142，Δp = 0.3536，Δx·Δp = 0.5000（高斯取等）")
    func minUncertainty() {
        let (x, psi2, p, phi2) = Self.state()
        let (_, dxStd) = OperatorRepresentationMath.moments(grid: x, density: psi2)
        let (_, dpStd) = OperatorRepresentationMath.moments(grid: p, density: phi2)
        #expect(abs(dxStd - 1.4142) < 1e-3, "Δx = \(dxStd)")
        #expect(abs(dpStd - 0.3536) < 1e-3, "Δp = \(dpStd)")
        #expect(abs(dxStd * dpStd - 0.5) < 1e-3, "乘积 = \(dxStd * dpStd)")
    }

    @Test("闭式对拍：|ψ(x)|² = (σ√π)⁻¹e^{−x²/σ²}（rel < 1e-9）")
    func positionClosedForm() {
        let (x, psi2, _, _) = Self.state()
        let s = Self.sigmaX
        for i in x.indices {
            let closed = exp(-x[i] * x[i] / (s * s)) / (s * .pi.squareRoot())
            #expect(abs(psi2[i] - closed) / closed < 1e-9,
                    "第\(i)点 \(psi2[i]) vs \(closed)")
        }
    }

    @Test("闭式对拍：|φ(p)|² = (σ/√π)e^{−σ²p²}，信号区 rel < 1e-6；|p| > 3 噪声底 < 1e-16")
    func momentumClosedForm() {
        let (_, _, p, phi2) = Self.state()
        let s = Self.sigmaX
        var checked = 0
        for i in p.indices {
            let ap = abs(p[i])
            if ap <= 1.5 {
                let closed = s / .pi.squareRoot() * exp(-s * s * p[i] * p[i])
                #expect(abs(phi2[i] - closed) / closed < 1e-6,
                        "第\(i)点 p=\(p[i]) \(phi2[i]) vs \(closed)")
                checked += 1
            }
        }
        #expect(checked >= 19, "信号区样本数 \(checked)")
        for i in p.indices where abs(p[i]) > 3.0 {
            #expect(phi2[i] < 1e-16, "噪声区 p=\(p[i]) y=\(phi2[i])")
        }
    }

    @Test("fixture 对拍：位置曲线 512 点最近网格（rel < 1e-7）")
    func positionMatchesFixture() throws {
        let fx = try ModuleFixture.load("算符表象_高斯波包傅里叶对偶__v2022")
        let (fxp, fyp) = try fx.line(0, 0, index: 0)
        let (x, psi2, _, _) = Self.state()
        let dx = x[1] - x[0]
        #expect(fxp.count == 512, "fixture 位置曲线点数 \(fxp.count)")
        for i in fxp.indices {
            let j = Int(((fxp[i] - x[0]) / dx).rounded())
            #expect(j >= 0 && j < x.count, "第\(i)点 x=\(fxp[i]) 越界")
            guard j >= 0 && j < x.count else { continue }
            #expect(abs(psi2[j] - fyp[i]) / max(fyp[i], 1e-300) < 1e-7,
                    "第\(i)点 \(psi2[j]) vs \(fyp[i])")
        }
    }

    @Test("fixture 对拍：动量信号点（y > 1e-16，|x| 反推网格，rel < 1e-7）")
    func momentumMatchesFixture() throws {
        let fx = try ModuleFixture.load("算符表象_高斯波包傅里叶对偶__v2022")
        let (fxp, fyp) = try fx.line(0, 1, index: 0)
        let (_, _, p, phi2) = Self.state()
        let n = p.count
        let dp = p[n / 2 + 1]        // p[N/2] = 0，p[N/2+1] = Δp
        var strong = 0, noise = 0
        for i in fxp.indices {
            guard fyp[i] > 1e-16 else { noise += 1; continue }
            let j = n / 2 + Int((abs(fxp[i]) / dp).rounded())
            #expect(j >= 0 && j < n, "第\(i)点 x=\(fxp[i]) 越界")
            guard j >= 0 && j < n else { continue }
            #expect(abs(phi2[j] - fyp[i]) / fyp[i] < 1e-7,
                    "信号点 第\(i)点 \(phi2[j]) vs \(fyp[i])")
            strong += 1
        }
        #expect(strong >= 4, "信号点数 \(strong)")
        #expect(strong + noise == fxp.count)
    }

    @Test("compute 输出结构：双图（512 点位置 + 截窗动量）+ 摘要 4 项 + 理论卡")
    func computeStructure() async throws {
        let module = OperatorRepresentationModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        guard case .lineSeries(let pos) = result.charts[0],
              case .lineSeries(let mom) = result.charts[1] else {
            Issue.record("应为两个 lineSeries"); return
        }
        #expect(pos.series.count == 1)
        #expect(pos.series[0].points.count == 512)      // 4096 / 8
        #expect(mom.series.count == 1)
        #expect(mom.series[0].points.count >= 32)       // |p| ≤ 4 截窗
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("秒级档预算：compute 中位 < 2000 ms")
    func computeBudget() async throws {
        let module = OperatorRepresentationModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 2000)
    }
}
