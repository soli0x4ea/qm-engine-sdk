import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W4-1 量子行为三实验对照：复刻脚本内 5 条断言（fixtures 量子行为_三实验对照__v2022.json）。
@Suite("W4 量子行为三实验对照")
struct QuantumBehaviorModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    /// 脚本默认几何：d=1, a=0.4, x ∈ [-3,3] × 3000。
    private func distributions() -> (x: [Double], bullets: [Double], waveB: [Double],
                                     electron: [Double], watched: [Double]) {
        QuantumBehaviorMath.distributions(d: 1.0, a: 0.4)
    }

    @Test("check1：电子主极大在 x=0，一级极小在 x=+0.5（条纹间距之半）")
    func electronExtrema() throws {
        let (x, _, _, electron, _) = distributions()
        // 主极大
        let iMax = electron.indices.max(by: { electron[$0] < electron[$1] })!
        #expect(abs(x[iMax]) < 0.01, "主极大 x=\(x[iMax])")
        // 一级极小（x>0 侧局部极小）
        var minsRight: [Double] = []
        for i in 1..<(x.count - 1) where electron[i] < electron[i - 1] && electron[i] < electron[i + 1] {
            if x[i] > 0 { minsRight.append(x[i]) }
        }
        let firstMin = try #require(minsRight.min())
        #expect(abs(firstMin - 0.5) < 0.05, "一级极小 x=\(firstMin)")
    }

    @Test("check2：被监视情形中央区（|x|<2）无归零深谷（min/peak > 0.05）")
    func watchedNoDeepDip() {
        let (x, _, _, _, watched) = distributions()
        let peak = watched.max()!
        let central = zip(x, watched).filter { abs($0.0) < 2.0 }
        let frac = central.map(\.1).min()! / peak
        #expect(frac > 0.05, "中央区 min/peak=\(frac)")
    }

    @Test("check3：电子极小深谷 < 1e-3，同位置监视情形 > 0.05")
    func dipContrast() throws {
        let (x, _, _, electron, watched) = distributions()
        // 电子一级极小位置（x>0）
        var minsRight: [Double] = []
        for i in 1..<(x.count - 1) where electron[i] < electron[i - 1] && electron[i] < electron[i + 1] {
            if x[i] > 0 { minsRight.append(x[i]) }
        }
        let firstMin = try #require(minsRight.min())
        let iMin = x.indices.min(by: { abs(x[$0] - firstMin) < abs(x[$1] - firstMin) })!
        #expect(electron[iMin] < 1e-3, "电子深谷 \(electron[iMin])")
        #expect(watched[iMin] > 0.05, "监视同位置 \(watched[iMin])")
    }

    @Test("归一化 sinc 与物理常数：sinc(0)=1；50 kV 电子 λ ≈ 5.49 pm")
    func sincAndDeBroglie() throws {
        #expect(QuantumBehaviorMath.sinc(0) == 1)
        #expect(abs(QuantumBehaviorMath.sinc(1)) < 1e-15)   // sin(π)/π = 0
        let cs = try constants()
        let lam = try cs.value("h") / (2 * cs.value("m_e") * cs.value("e") * 50e3).squareRoot()
        // 非相对论 50 kV 电子德布罗意波长公认值 ≈ 5.4859 pm（脚本打印值 5.49 pm 量级）
        #expect(abs(lam * 1e12 - 5.4859) < 0.01, "λ=\(lam * 1e12) pm")
    }

    @Test("compute 输出结构：4 图 4 系列 + 参考线 + 摘要 + 理论卡")
    func computeStructure() async throws {
        let module = QuantumBehaviorModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 4)
        for c in result.charts {
            guard case .lineSeries(let l) = c else { Issue.record("应全为 lineSeries"); return }
            #expect(l.series.count == 1)
            #expect(l.series[0].points.count == 500)   // 3000 / stride 6
        }
        guard case .lineSeries(let electronChart) = result.charts[2] else { return }
        #expect(electronChart.referenceLines.count == 1)
        #expect(abs(electronChart.referenceLines[0].value - 0.5) < 1e-12)
        #expect(result.summary.count == 3)
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("实时档预算：compute 中位 < 16 ms")
    func computeBudget() async throws {
        let module = QuantumBehaviorModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }
}
