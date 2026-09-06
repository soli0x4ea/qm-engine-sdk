import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W5-10 BB84 密钥率：R(Q) = 1 − 2h(Q)（脚本 [check] 打印值断言 + fixtures 曲线对拍）。
@Suite("W5 量子密钥分发：BB84 密钥率")
struct BB84ModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    @Test("二元熵端点：h(0)=0、h(0.5)=1、h(1)=0；h(0.11) ≈ 0.49992（脚本打印值）")
    func entropyEndpoints() {
        #expect(BB84Math.binaryEntropy(0) == 0)
        #expect(BB84Math.binaryEntropy(1) == 0)
        #expect(abs(BB84Math.binaryEntropy(0.5) - 1) < 1e-12)
        // 脚本 [check] h(0.11)=0.49995 → 精确 0.499916（5 位打印舍入到 0.49992）
        #expect(abs(BB84Math.binaryEntropy(0.11) - 0.499916) < 1e-5, "h(0.11) = \(BB84Math.binaryEntropy(0.11))")
    }

    @Test("密钥率：R(0)=1、R(0.5%)≈0.909（脚本打印值）；阈值 11% 处归零")
    func keyRateValues() {
        #expect(abs(BB84Math.keyRate(0) - 1) < 1e-12)
        #expect(abs(BB84Math.keyRate(0.005) - 0.909165) < 1e-4, "R(0.5%) = \(BB84Math.keyRate(0.005))")
        // R(11%) = 1 − 2×0.499916 ≈ 1.68e-4：阈值处「~0」（脚本 [check] 口径）
        #expect(abs(BB84Math.keyRate(0.11)) < 1e-3)
        #expect(BB84Math.keyRate(0.11) > 0, "阈值处仍非负")
        #expect(BB84Math.keyRate(0.12) < 0, "12% 已不安全")
        // 单调下降 + 对称性 h(p)=h(1−p)
        #expect(BB84Math.keyRate(0.05) > BB84Math.keyRate(0.10))
        #expect(abs(BB84Math.binaryEntropy(0.2) - BB84Math.binaryEntropy(0.8)) < 1e-12)
    }

    @Test("fixture 曲线对拍：R(Q) 主曲线（512 点降采样，Q ∈ [0, 25]%）")
    func curveMatchFixture() throws {
        let fx = try ModuleFixture.load("量子密钥分发_BB84密钥率__v2022")
        let (x, y) = try fx.line(0, 0, index: 0)
        #expect(x.count == 512)
        #expect(abs(x.first! - 0) < 1e-12 && abs(x.last! - 25) < 1e-9)
        // R 在 Q≈11% 处 ~1.7e-4（小值），x 9 位舍入经 |R′(Q)|≈4 传导后 rel 放大 → 容差 5e-6
        expectPointwiseClose(x.map { BB84Math.keyRate($0 / 100) }, y,
                             tolerance: 5e-6, "R(Q) = 1 − 2h(Q)")
    }

    @Test("compute 输出结构：300 点主曲线 + 安全区填充（264 点 ≤ 11%）+ 2 参考线 + 摘要 + 协议示意图")
    func computeStructure() async throws {
        let module = BB84Module()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)   // W13d：+协议泳道示意图（笔记 32）
        guard case .lineSeries(let c0) = result.charts[0] else {
            Issue.record("应为 lineSeries"); return
        }
        #expect(c0.series.count == 1)
        #expect(c0.series[0].points.count == 300)   // 600 / 2
        #expect(c0.referenceLines.count == 2)
        #expect(c0.referenceLines.contains { $0.value == 11 && $0.axis == .x })
        // 安全区填充：Q ≤ 11% 共 264 点（步长 0.25/599），基线 0
        #expect(c0.areaFills.count == 1)
        let fill = c0.areaFills[0]
        #expect(fill.points.count == 264, "实际 \(fill.points.count)")
        #expect(fill.baseline == 0)
        #expect(fill.points.allSatisfy { $0.x <= 11 })
        #expect(abs(fill.points.last!.x - 10.976) < 1e-2)
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 4)
        // 协议示意图：4 泳道 + 12 过程标记点（4 Alice + 4 光子 + 4 Bob）+ 3 密钥点 + 箭头
        guard case .schematic(let s1) = result.charts[1] else {
            Issue.record("应为 schematic"); return
        }
        #expect(s1.bands.count == 4)
        #expect(s1.markers.count == 15)
        #expect(s1.markers.filter { !$0.filled }.count == 1, "轮 3 基矢误选 = 空心点")
        #expect(s1.callouts.count == 9)   // 8 泳道间箭头 + 1 基矢比对标注
    }

    @Test("实时档预算：compute 中位 < 16 ms")
    func computeBudget() async throws {
        let module = BB84Module()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }
}
