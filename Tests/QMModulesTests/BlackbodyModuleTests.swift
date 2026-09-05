import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W3-2 黑体三律模块：fixtures 对拍（黑体辐射_三律对比__v2022.json，T = 5772 K）。
///
/// fixture 曲线为 512 点重采样，插值误差实测 ~9×10⁻⁸，对拍容差取 1×10⁻⁶。
@Suite("W3 黑体三律模块")
struct BlackbodyModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    @Test("三律谱曲线逐点对拍 fixtures")
    func spectraMatchFixtures() throws {
        let fx = try ModuleFixture.load("黑体辐射_三律对比__v2022")
        let cs = try constants()
        let h = try cs.value("h"), kB = try cs.value("kB"), c = try cs.value("c")
        let T = 5772.0

        let (px, py) = try fx.line(0, 0, label: "Planck")
        expectPointwiseClose(
            px.map { BlackbodyMath.uLambda($0 * 1e-9, T, h: h, kB: kB, c: c) * 1e9 },
            py, tolerance: 1e-6, "普朗克谱")

        let (wx, wy) = try fx.line(0, 0, label: "Wien")
        expectPointwiseClose(
            wx.map { BlackbodyMath.uWien($0 * 1e-9, T, h: h, kB: kB, c: c) * 1e9 },
            wy, tolerance: 1e-6, "维恩近似谱")

        let (rx, ry) = try fx.line(0, 0, label: "Rayleigh")
        expectPointwiseClose(
            rx.map { BlackbodyMath.uRayleighJeans($0 * 1e-9, T, kB: kB, c: c) * 1e9 },
            ry, tolerance: 1e-6, "瑞利-金斯谱")

        // RJ 截断规则（u_RJ < 4×普朗克峰值）起点与 fixture 一致
        let grid = Num.linspace(50, 3000, count: 4000)
        let peak = grid.map {
            BlackbodyMath.uLambda($0 * 1e-9, T, h: h, kB: kB, c: c) * 1e9
        }.max() ?? 0
        let firstMasked = grid.first {
            BlackbodyMath.uRayleighJeans($0 * 1e-9, T, kB: kB, c: c) * 1e9 < 4 * peak
        }
        #expect(firstMasked.map { abs($0 - rx[0]) < 0.01 } == true,
                "RJ 截断起点 \(String(describing: firstMasked)) vs fixture \(rx[0])")
    }

    @Test("斯特藩-玻尔兹曼与维恩位移数值验证（复刻脚本内断言）")
    func stefanAndWienChecks() throws {
        let cs = try constants()
        let h = try cs.value("h"), kB = try cs.value("kB"), c = try cs.value("c")
        let T = 5772.0

        // check1：全频段积分 j* = c/4 · ∫u_ν dν vs σT⁴（脚本断言 rel < 1e-3）
        let nus = Num.linspace(1e9, 60 * kB * T / h, count: 60000)
        let U = Num.trapezoid(nus, nus.map { BlackbodyMath.uNu($0, T, h: h, kB: kB, c: c) })
        let sigma = BlackbodyMath.sigmaDerived(h: h, kB: kB, c: c)
        let jAna = sigma * pow(T, 4)
        let rel1 = abs(c / 4 * U - jAna) / jAna
        #expect(rel1 < 1e-3, "斯特藩-玻尔兹曼积分 rel=\(rel1)")
        // σ 导出值与 CODATA 2022 推荐值（5.670374419e-8）自洽
        #expect(abs(sigma - 5.670374419e-8) / 5.670374419e-8 < 1e-6)

        // check2：λmax 网格搜索 vs b/T（脚本断言 rel < 1e-3）
        let lams = Num.linspace(1e-8, 8e-6, count: 200000)
        var peakVal = 0.0, peakLam = 0.0
        for lam in lams {
            let u = BlackbodyMath.uLambda(lam, T, h: h, kB: kB, c: c)
            if u > peakVal { peakVal = u; peakLam = lam }
        }
        let b = BlackbodyMath.wienB(h: h, kB: kB, c: c)
        let rel2 = abs(peakLam - b / T) / (b / T)
        #expect(rel2 < 1e-3, "维恩位移峰值 rel=\(rel2)")
    }

    @Test("compute 输出结构：系列 / 参考线 / 摘要 / 理论卡")
    func computeStructure() async throws {
        let module = BlackbodyModule()
        let values = ParamValues.defaults(for: module.params)
        let result = try await module.compute(values, constants: try constants())

        guard case .lineSeries(let line) = result.charts[0] else {
            Issue.record("charts[0] 应为 lineSeries"); return
        }
        #expect(line.series.map(\.name) == ["普朗克", "维恩近似", "瑞利-金斯"])
        // 普朗克/维恩：4000 点网格 stride 8 → 500 点；
        // 瑞利-金斯：截断（u_RJ < 4×峰值）后 ~2951 点再抽稀 → ~369 点
        #expect(line.series[0].points.count == 500)
        #expect(line.series[1].points.count == 500)
        #expect((100..<500).contains(line.series[2].points.count),
                "RJ 系列点数 \(line.series[2].points.count)")
        #expect(line.referenceLines.count == 1)

        let cs = try constants()
        let lamMax = BlackbodyMath.wienB(h: try cs.value("h"), kB: try cs.value("kB"),
                                         c: try cs.value("c")) / 5772 * 1e9
        #expect(abs(line.referenceLines[0].value - lamMax) < 1e-9, "λmax 参考线")
        #expect(result.summary.count == 3)
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("实时档预算：compute 中位 < 16 ms")
    func computeBudget() async throws {
        let module = BlackbodyModule()
        try await expectComputeUnderBudget(
            module: module,
            values: ParamValues.defaults(for: module.params),
            constants: try constants(),
            budgetMillis: 16)
    }
}
