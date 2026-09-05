import Testing
import Foundation
import Numerics
@testable import QMModules
@testable import EngineKit

/// W6 笔记 34《量子计算与量子算法》：QFT 8×8 酉矩阵（swift-numerics 复数通道）。
/// fixtures：量子计算_QFT__v2022.json（输入 |1⟩ 的 8 柱输出概率，每柱 1/8）。
@Suite("W6 量子傅里叶变换")
struct QFTModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    // MARK: - 物理律（酉性 / 谱性质）

    @Test("物理律：酉性 U†U = I（n = 2/3/4，max 偏差 < 1e-12）")
    func unitarityLaw() {
        for n in [2, 3, 4] {
            let err = QFTMath.unitaryError(QFTMath.qftMatrix(n: n))
            #expect(err < 1e-12, "n=\(n)：max|U†U − I| = \(err)")
        }
        // 门级酉性：受控相位门与嵌入
        #expect(QFTMath.unitaryError(QFTMath.controlledPhase(k: 2)) < 1e-12)
        #expect(QFTMath.unitaryError(QFTMath.controlledPhase(k: 3)) < 1e-12)
        #expect(QFTMath.unitaryError(
            QFTMath.embed2q(QFTMath.controlledPhase(k: 2), n: 3, control: 0, target: 1)) < 1e-12)
        #expect(QFTMath.unitaryError(QFTMath.hadamardEmbedded(n: 3, q: 1)) < 1e-12)
    }

    @Test("物理律：QFT² = 比特取反置换（谱性质，< 1e-12）；元素等幅 |U[y][x]|² = 1/N")
    func spectralLaw() {
        for n in [2, 3, 4] {
            #expect(QFTMath.squareError(n: n) < 1e-12, "n=\(n)：QFT² = P_neg")
            #expect(QFTMath.columnNormError(n: n) < 1e-14, "n=\(n)：元素模方恒 1/N")
        }
        // 取反置换自检：P_neg² = I（取反是 对合）
        let p = QFTMath.negationPermutation(n: 3)
        #expect(QFTMath.unitaryError(QFTMath.matmul(p, p)) < 1e-15)
    }

    @Test("物理律：QFT|0⟩ = 均匀叠加（等幅 1/√8、零相位）；定义式逐元素")
    func definitionLaw() throws {
        let psi0 = QFTMath.outputAmplitudes(n: 3, inputIndex: 0)
        #expect(psi0.count == 8)
        for z in psi0 {
            #expect(abs(z.real - 1.0 / 8.0.squareRoot()) < 1e-15)
            #expect(abs(z.imaginary) < 1e-15)
        }
        // 输出概率归一
        let probs = QFTMath.outputProbabilities(n: 3, inputIndex: 1)
        #expect(abs(probs.reduce(0, +) - 1) < 1e-12)

        // 定义式 U[y][x] = e^{2πi·xy/8}/√8 的代表元素：
        // U[1][1] = e^{iπ/4}/√8；U[4][1] = e^{iπ}/√8 = −1/√8；U[2][3] = e^{3πi/2}/√8 = −i/√8
        let u = QFTMath.qftMatrix(n: 3)
        let norm = 1.0 / 8.0.squareRoot()
        #expect(abs(u[1][1].real - norm * 2.0.squareRoot() / 2.0) < 1e-15)
        #expect(abs(u[1][1].imaginary - norm * 2.0.squareRoot() / 2.0) < 1e-15)
        #expect(abs(u[4][1].real + norm) < 1e-15 && abs(u[4][1].imaginary) < 1e-15)
        #expect(abs(u[2][3].imaginary + norm) < 1e-15 && abs(u[2][3].real) < 1e-15)
    }

    // MARK: - fixtures 对拍

    @Test("fixtures：输入 |1⟩ 的 8 柱输出概率逐柱 1/8（< 1e-8）")
    func barsMatchFixtures() throws {
        let fx = try ModuleFixture.load("量子计算_QFT__v2022")
        let bars = try fx.bars(0, 0)
        #expect(bars.count == 8)
        let probs = QFTMath.outputProbabilities(n: 3, inputIndex: 1)
        for i in 0..<8 {
            #expect(bars[i].x == Double(i))
            #expect(abs(bars[i].height - 1.0 / 8.0) < 1e-8,
                    "柱 \(i)：\(bars[i].height) vs 1/8")
            #expect(abs(bars[i].height - probs[i]) < 1e-12,
                    "与 Swift 概率一致（等幅是傅里叶变换的标志）")
        }
    }

    // MARK: - compute 结构与预算

    @Test("compute 结构：柱图 + Re/Im 双曲线 + 酉性摘要五卡；秒级档预算 < 2 s")
    func computeStructureAndBudget() async throws {
        let module = QFTModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        guard case .bars(let bars) = result.charts[0] else {
            Issue.record("charts[0] 应为 bars"); return
        }
        #expect(bars.bars.count == 8)
        for b in bars.bars { #expect(abs(b.value - 0.125) < 1e-12) }
        guard case .lineSeries(let reIm) = result.charts[1] else {
            Issue.record("charts[1] 应为 lineSeries"); return
        }
        #expect(reIm.series.count == 2)
        #expect(reIm.series[0].points.count == 8 && reIm.series[1].points.count == 8)
        #expect(result.summary.count == 5)
        #expect(result.theory?.formulas.count == 4)

        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 2000)
    }
}
