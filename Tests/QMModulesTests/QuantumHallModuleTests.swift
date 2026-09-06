import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W5-8 量子霍尔：Landau 能级 + IQHE 平台（脚本打印值断言）
/// + fixtures 双 Y 轴两曲线对拍。
@Suite("W5 量子霍尔：Landau 能级与平台")
struct QuantumHallModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    /// 脚本默认：n_e = 2.4e15 m⁻²（二维电子气典型值）。
    private let ne = 2.4e15

    @Test("R_K = h/e² = 25812.8 Ω（脚本打印值）")
    func vonKlitzing() throws {
        let cs = try constants()
        let h = try cs.value("h"), e = try cs.value("e")
        let rk = h / (e * e)
        #expect(abs(rk / 25812.8 - 1) < 1e-6, "R_K = \(rk) Ω")
        #expect(abs(rk / 25812.80745 - 1) < 1e-7, "2019 SI 精确值")
    }

    @Test("Landau 能级：ℏω_c @10 T = 1.158 meV；简并度 eB/h = 2.418e15 m⁻²")
    func landau() throws {
        let cs = try constants()
        let hbar = try cs.value("hbar"), me = try cs.value("m_e"),
            e = try cs.value("e"), eV = try cs.value("eV"), h = try cs.value("h")
        // E_0 = ½ℏω_c（脚本口径：ℏω_c = ℏeB/m_e）
        let wc10 = 2.0 * QuantumHallMath.landauLevel(10, n: 0, hbar: hbar, me: me, e: e) / eV * 1e3
        #expect(abs(wc10 - 1.158) < 1e-3, "ℏω_c = \(wc10) meV")
        #expect(abs(e * 10 / h / 2.418e15 - 1) < 1e-3, "eB/h")
        // 等间距：E_n − E_{n−1} = ℏω_c（脚本口径，非 E_0）
        let spacing = QuantumHallMath.landauLevel(10, n: 1, hbar: hbar, me: me, e: e)
            - QuantumHallMath.landauLevel(10, n: 0, hbar: hbar, me: me, e: e)
        #expect(abs(spacing / (hbar * e * 10 / me) - 1) < 1e-12)
        for n in 1...5 {
            let d = QuantumHallMath.landauLevel(10, n: n, hbar: hbar, me: me, e: e)
                - QuantumHallMath.landauLevel(10, n: n - 1, hbar: hbar, me: me, e: e)
            #expect(abs(d - spacing) < 1e-30)
        }
    }

    @Test("IQHE 玩具模型：ν 与 ρ_xy 平台（脚本打印四档）")
    func plateau() throws {
        let cs = try constants()
        let h = try cs.value("h"), e = try cs.value("e")
        // 脚本四档：Nfilled = max(1, ⌊ν⌋)（B=10 T 时 ν≈0.99 被钳到 1）
        let expected: [(Double, Double, Double)] = [
            (10.0, 1, 25812.8), (5.0, 1, 25812.8), (3.33, 2, 12906.4), (2.0, 4, 6453.2),
        ]
        for (b, nFill, rExp) in expected {
            let nu = QuantumHallMath.fillingFactor(b, ne: ne, h: h, e: e)
            #expect(abs(QuantumHallMath.filledLevels(nu) - Double(nFill)) < 0.01,
                   "B=\(b) T ν=\(nu) N=\(QuantumHallMath.filledLevels(nu))")
            let r = QuantumHallMath.rhoXY(b, ne: ne, h: h, e: e)
            #expect(abs(r / rExp - 1) < 1e-3, "B=\(b) T ρ_xy = \(r)")
        }
        // 精确整数填充 B = ν·(h/(e·ν)) 处 dist = 0 → ρ_xx = 0（平台中心）
        let bNu1 = ne * h / (e * 1.0)
        let bNu2 = ne * h / (e * 2.0)
        #expect(QuantumHallMath.rhoXX(bNu1, ne: ne, h: h, e: e) < 1e-9)
        #expect(QuantumHallMath.rhoXX(bNu2, ne: ne, h: h, e: e) < 1e-9)
        // 精确半填充 dist = 0.5 → 0.45 截断饱和值 1000
        let bHalf = ne * h / (e * 1.5)
        #expect(abs(QuantumHallMath.rhoXX(bHalf, ne: ne, h: h, e: e) - 1000) < 1e-6)
    }

    @Test("fixture 曲线对拍：6 条 Landau 能级 + ρ_xy/R_K + ρ_xx")
    func curvesMatchFixtures() throws {
        let fx = try ModuleFixture.load("量子霍尔_Landau能级与平台__v2022")
        let cs = try constants()
        let hbar = try cs.value("hbar"), me = try cs.value("m_e"),
            e = try cs.value("e"), eV = try cs.value("eV"), h = try cs.value("h")
        let rk = h / (e * e)

        // fig0：E_n(B)（B = linspace(0, 15, 400) 的降采样）
        for n in 0..<6 {
            let (x, y) = try fx.line(0, 0, index: n)
            #expect(x.count == 400)
            expectPointwiseClose(
                x.map { QuantumHallMath.landauLevel($0, n: n, hbar: hbar, me: me, e: e) / eV * 1e3 },
                y, tolerance: 1e-8, "E_\(n)(B)")
        }
        // fig1 ax0：ρ_xy/R_K（B2 = linspace(1, 10, 1200) 的降采样）
        let (x1, y1) = try fx.line(1, 0, index: 0)
        #expect(x1.count == 512)
        expectPointwiseClose(x1.map { QuantumHallMath.rhoXY($0, ne: ne, h: h, e: e) / rk },
                             y1, tolerance: 1e-8, "ρ_xy/R_K")
        // fig1 ax1：ρ_xx（分段线性 + 0.45 截断拐点导数跳变，x 9 位舍入在拐点
        // 附近放大两个量级 → 容差 1e-5；W4 先例：导数陡峭处按函数特性分层）
        let (x2, y2) = try fx.line(1, 1, index: 0)
        #expect(x2.count == 512)
        expectPointwiseClose(x2.map { QuantumHallMath.rhoXX($0, ne: ne, h: h, e: e) },
                             y2, tolerance: 1e-5, "ρ_xx (dips)")
    }

    @Test("compute 输出结构：Landau 6 系列 + 双 Y 轴 1+1 系列")
    func computeStructure() async throws {
        let module = QuantumHallModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        guard case .lineSeries(let c0) = result.charts[0],
              case .dualAxisLineSeries(let c1) = result.charts[1] else {
            Issue.record("应为 lineSeries + dualAxisLineSeries"); return
        }
        #expect(c0.series.count == 6)
        #expect(c0.series.allSatisfy { $0.points.count == 400 })
        #expect(c1.primary.count == 1)
        #expect(c1.secondary.count == 1)
        #expect(c1.primary[0].points.count == 600)   // 1200 / 2
        #expect(c1.secondary[0].points.count == 600)
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("实时档预算：compute 中位 < 16 ms")
    func computeBudget() async throws {
        let module = QuantumHallModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }
}
