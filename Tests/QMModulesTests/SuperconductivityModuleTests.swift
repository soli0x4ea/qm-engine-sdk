import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W5-7 超导：BCS 能隙插值 + Fraunhofer 调制（脚本打印值断言；
/// 两主曲线无 label 未入 fixtures，公式复算对拍——W4 双缝同款口径）。
@Suite("W5 超导：BCS 与约瑟夫森")
struct SuperconductivityModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    @Test("能隙插值：Δ(0.5T_c)/Δ(0) = tanh(1.76) = 0.9425（脚本打印 0.9428 为网格点 0.4992 处值）")
    func gapHalf() {
        #expect(abs(SuperconductivityMath.gapRatio(0.5) - tanh(1.76)) < 1e-12)
        #expect(abs(SuperconductivityMath.gapRatio(0.5) - 0.9425) < 1e-4)
        // 脚本网格 Tr = linspace(0.02, 1.4, 600) 最近 0.5 的点为 0.49919 → 打印 0.9428
        #expect(abs(SuperconductivityMath.gapRatio(0.49919) - 0.9428) < 1e-3)
    }

    @Test("能隙端点：T→0 平台 →1；T ≥ T_c 闭合为 0")
    func gapEndpoints() {
        #expect(abs(SuperconductivityMath.gapRatio(0.001) - 1) < 1e-3)
        #expect(SuperconductivityMath.gapRatio(1) == 0)
        #expect(SuperconductivityMath.gapRatio(1.4) == 0)
        // T_c 附近 √ 收口：gapRatio(0.99) 微小但非零
        let g99 = SuperconductivityMath.gapRatio(0.99)
        #expect(g99 > 0 && g99 < 0.25)
    }

    @Test("Δ(0) = 1.76 k_BT_c：Nb（T_c = 9.25 K）≈ 1.403 meV")
    func delta0() throws {
        let cs = try constants()
        let kB = try cs.value("kB"), eV = try cs.value("eV")
        let d0 = 1.76 * kB * 9.25
        #expect(abs(d0 / eV * 1e3 - 1.4029) < 1e-3, "Δ(0) = \(d0 / eV * 1e3) meV")
    }

    @Test("Fraunhofer：中心 1、整数磁通过零、|Φ/Φ₀| ≤ 0.5 内 sinc 包络")
    func fraunhofer() {
        #expect(abs(SuperconductivityMath.fraunhofer(0) - 1) < 1e-12)
        for k in [-3, -2, -1, 1, 2, 3] {
            #expect(abs(SuperconductivityMath.fraunhofer(Double(k))) < 1e-12, "Φ = \(k)Φ₀")
        }
        #expect(abs(SuperconductivityMath.fraunhofer(0.5) - 2 / .pi) < 1e-12)
        // 包络：任意 Φ 处 |sinc| ≤ 1/(πΦ/Φ₀)
        for x in stride(from: -2.9, through: 2.9, by: 0.37) {
            #expect(SuperconductivityMath.fraunhofer(x) <= 1 / (.pi * abs(x)) + 1e-12)
        }
        // 对称：I_c(−Φ) = I_c(+Φ)
        #expect(abs(SuperconductivityMath.fraunhofer(0.83)
                    - SuperconductivityMath.fraunhofer(-0.83)) < 1e-15)
    }

    @Test("compute 输出结构：2 图（300 + 450 点）+ 参考线 + 摘要 + 理论卡")
    func computeStructure() async throws {
        let module = SuperconductivityModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        guard case .lineSeries(let c0) = result.charts[0],
              case .lineSeries(let c1) = result.charts[1] else {
            Issue.record("应为 lineSeries"); return
        }
        #expect(c0.series[0].points.count == 300)   // 600 / 2
        #expect(c1.series[0].points.count == 450)   // 900 / 2
        #expect(c0.referenceLines.count == 1)
        #expect(c1.referenceLines.count == 2)
        #expect(result.summary.count == 3)
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("实时档预算：compute 中位 < 16 ms")
    func computeBudget() async throws {
        let module = SuperconductivityModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }
}
