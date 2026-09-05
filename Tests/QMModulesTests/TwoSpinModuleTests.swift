import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W4-16 两自旋单态三重态：复刻脚本断言（ρA = 𝟙/2、S = ln 2）+ fixtures 双线对拍。
/// 热图（imshow）未被导出管线捕获，数值由闭式断言覆盖。
@Suite("W4 两自旋单态三重态")
struct TwoSpinModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    private static func rhoA(_ stateID: String) -> [[TwoSpinMath.Z]] {
        TwoSpinMath.partialTraceOverB(TwoSpinMath.outer(TwoSpinMath.bellState(stateID)))
    }

    @Test("单态约化密度矩阵 ρA = 𝟙/2（脚本打印矩阵，|偏差| < 1e-15）")
    func singletReduced() {
        let r = Self.rhoA("S")
        #expect(abs(r[0][0].re - 0.5) < 1e-15 && abs(r[0][0].im) < 1e-15)
        #expect(abs(r[1][1].re - 0.5) < 1e-15 && abs(r[1][1].im) < 1e-15)
        #expect(r[0][1].abs2 < 1e-30 && r[1][0].abs2 < 1e-30)
    }

    @Test("四态熵：S/T0 = ln 2（1 ebit）；T± = 0（直积态）")
    func entropies() {
        let sS = TwoSpinMath.vonNeumannEntropy(Self.rhoA("S"))
        let sT0 = TwoSpinMath.vonNeumannEntropy(Self.rhoA("T0"))
        let sTp = TwoSpinMath.vonNeumannEntropy(Self.rhoA("Tp"))
        let sTm = TwoSpinMath.vonNeumannEntropy(Self.rhoA("Tm"))
        #expect(abs(sS - log(2.0)) < 1e-12, "S(单态) = \(sS)")
        #expect(abs(sT0 - log(2.0)) < 1e-12, "S(T0) = \(sT0)")
        #expect(sTp < 1e-15, "S(T+) = \(sTp)")
        #expect(sTm < 1e-15, "S(T−) = \(sTm)")
        // 脚本打印行：entropy ln2 nat = 1.0000 ebit
        #expect(abs(sS / log(2.0) - 1.0) < 1e-12)
    }

    @Test("T± 直积投影：ρA = |↑⟩⟨↑| / |↓⟩⟨↓|")
    func productStates() {
        let rp = Self.rhoA("Tp"), rm = Self.rhoA("Tm")
        #expect(abs(rp[0][0].re - 1) < 1e-15 && rp[1][1].abs2 < 1e-30)
        #expect(abs(rm[1][1].re - 1) < 1e-15 && rm[0][0].abs2 < 1e-30)
        #expect(rp[0][1].abs2 < 1e-30 && rm[0][1].abs2 < 1e-30)
    }

    @Test("fixture 对拍：关联双线 200 点（∓cos θ，rel < 1e-8）")
    func correlationsMatchFixtures() throws {
        let fx = try ModuleFixture.load("角动量与自旋_两自旋单态三重态__v2022")
        let (sx, sy) = try fx.line(0, 1, label: "singlet")
        let (tx, ty) = try fx.line(0, 1, label: "triplet")
        #expect(sx.count == 200 && tx.count == 200)
        // 脚本网格 θᵢ = i·π/199（linspace(0, π, 200)）；fixture x 只存 9 位有效数字，
        // 在 cos 过零附近放大相对误差，故闭式按精确网格复算。
        for i in sx.indices {
            let theta = .pi * Double(i) / 199.0
            #expect(abs(sx[i] - theta) / .pi < 1e-8, "θ 网格 第\(i)点 \(sx[i])")
            #expect(abs(tx[i] - theta) / .pi < 1e-8, "θ 网格 第\(i)点 \(tx[i])")
            let closedS = -cos(theta)
            let closedT = cos(theta)
            #expect(abs(sy[i] - closedS) / max(abs(closedS), 1e-300) < 1e-8,
                    "singlet 第\(i)点 \(sy[i]) vs \(closedS)")
            #expect(abs(ty[i] - closedT) / max(abs(closedT), 1e-300) < 1e-8,
                    "triplet 第\(i)点 \(ty[i]) vs \(closedT)")
        }
        // 端点语义：θ = 0 反平行 −1 / 平行 +1；θ = π/2 两条线过零
        #expect(abs(sy[0] + 1) < 1e-9 && abs(ty[0] - 1) < 1e-9)
        #expect(abs(sy[99]) < 1e-2 && abs(ty[99]) < 1e-2)
    }

    @Test("compute 输出结构：热图（2×2）+ 双线图 + 摘要 4 项 + 理论卡")
    func computeStructure() async throws {
        let module = TwoSpinModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        guard case .heatmap(let hm) = result.charts[0],
              case .lineSeries(let line) = result.charts[1] else {
            Issue.record("应为 heatmap + lineSeries"); return
        }
        #expect(hm.values.count == 2 && hm.values[0].count == 2)
        #expect(hm.xTicks == ["|↑⟩", "|↓⟩"] && hm.yTicks.count == 2)
        // 默认单态：热图 = 𝟙/2
        #expect(abs(hm.values[0][0] - 0.5) < 1e-15)
        #expect(abs(hm.values[1][1] - 0.5) < 1e-15)
        #expect(hm.values[0][1] == 0 && hm.values[1][0] == 0)
        #expect(line.series.count == 2)
        #expect(line.series[0].points.count == 200)
        #expect(line.series[1].points.count == 200)
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("秒级档预算：compute 中位 < 2000 ms")
    func computeBudget() async throws {
        let module = TwoSpinModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 2000)
    }
}
