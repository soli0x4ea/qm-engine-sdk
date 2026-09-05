import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W6 笔记 40《量子场论的基本图像》：Casimir 压强 / 一维振子链模式量子化。
/// fixtures：量子场论的基本图像_Casimir__v2022.json（无标签主曲线被 runner 过滤，
/// 仅保留 ideal exponent −4 参考线——方法验证由物理律测试承担）、
/// 量子场论的基本图像_模式量子化__v2022.json（BE 占据双线各 50 点；前两图无标签被过滤）。
@Suite("W6 量子场论基本图像")
struct QFTBasicsModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    // MARK: - 物理律（Casimir）

    @Test("物理律：吸引压强 P < 0 且严格 ∝ d⁻⁴——P·d⁴ 恒等于 π²ℏc/240")
    func casimirPowerLaw() throws {
        let c = try constants()
        let hbar = try c.value("hbar")
        let lightSpeed = try c.value("c")
        let coeff = Double.pi * Double.pi * hbar * lightSpeed / 240.0

        for d in [1e-7, 3.7e-7, 1e-6, 5.3e-6, 1e-5] {
            let p = QFTBasicsMath.casimirPressure(d, hbar: hbar, c: lightSpeed)
            #expect(QFTBasicsMath.casimirPressureSigned(d, hbar: hbar, c: lightSpeed) < 0,
                    "吸引：带符号压强恒负")
            #expect(abs(p * pow(d, 4) - coeff) / coeff < 1e-12, "d = \(d)：P·d⁴ = 常数")
        }
        // 间距 ×10 ⇒ 压强 ÷10⁴（d⁻⁴ 幂律的直观形式）
        let p1 = QFTBasicsMath.casimirPressure(1e-6, hbar: hbar, c: lightSpeed)
        let p01 = QFTBasicsMath.casimirPressure(1e-7, hbar: hbar, c: lightSpeed)
        let p10 = QFTBasicsMath.casimirPressure(1e-5, hbar: hbar, c: lightSpeed)
        #expect(abs(p01 / p1 - 1e4) < 1e-8)
        #expect(abs(p1 / p10 - 1e4) < 1e-8)
    }

    @Test("物理律：已知量级——d = 1 μm 处 |F|/A ≈ 1.3 mPa；A = 1 cm² 吸引力 ≈ 1.3×10⁻⁷ N")
    func casimirMagnitudeLaw() throws {
        let c = try constants()
        let hbar = try c.value("hbar")
        let lightSpeed = try c.value("c")
        let p1um = QFTBasicsMath.casimirPressure(1e-6, hbar: hbar, c: lightSpeed)
        #expect(p1um > 1.29e-3 && p1um < 1.31e-3, "|F|/A = \(p1um) N/m² ≈ 1.3 mPa")
        let force1cm2 = p1um * 1e-4
        #expect(force1cm2 > 1.29e-7 && force1cm2 < 1.31e-7,
                "A = 1 cm²：F = \(force1cm2) N ≈ 1.3×10⁻⁷ N（约一个红细胞的重量）")
    }

    @Test("物理律：局部幂律指数 d ln P / d ln d ≡ −4（解析幂律的全局不变量）")
    func localExponentLaw() throws {
        let c = try constants()
        let hbar = try c.value("hbar")
        let lightSpeed = try c.value("c")
        // 模块同口径：0.1 ~ 10 μm 对数网格 200 点
        let d = Num.linspace(log10(1e-7), log10(1e-5), count: 200).map { pow(10, $0) }
        let p = d.map { QFTBasicsMath.casimirPressure($0, hbar: hbar, c: lightSpeed) }
        let slopes = QFTBasicsMath.localExponents(d, p)
        #expect(slopes.count == 199)
        for s in slopes { #expect(abs(s + 4) < 1e-12, "局部指数 = \(s)") }
    }

    // MARK: - 物理律（模式量子化）

    @Test("物理律：声学支 ω(k→0) → v|k|（线性色散）；布里渊区边界达带宽 2v/a")
    func acousticDispersionLaw() throws {
        _ = try constants()  // 预热常量缓存（本测试用解析值，不消费 CODATA 项）
        let a = 5e-10, v = 1e3
        let ks = QFTBasicsMath.waveVectors(N: 50, a: a)
        let omegas = ks.map { QFTBasicsMath.dispersion($0, a: a, v: v) }
        // k → 0 线性（最小非零 |k| 模：|sin(x)/x − 1| < 1e-3）
        let iSmall = ks.indices.filter { omegas[$0] > 0 }
            .min { abs(ks[$0]) < abs(ks[$1]) } ?? 0
        #expect(abs(omegas[iSmall] / abs(ks[iSmall]) - v) / v < 1e-3,
                "ω/k → v（声速）：\(omegas[iSmall] / abs(ks[iSmall])) vs \(v)")
        // 带宽上界：ω ≤ 2v/a，且布里渊区边界 m = ±25 处取等（|sin(π/2)| = 1）
        let omegaMax = 2.0 * v / a
        for (k, w) in zip(ks, omegas) {
            #expect(w <= omegaMax + 1e-6, "ω ≤ 2v/a（带宽）")
            if abs(abs(k) - .pi / a) < 1e-6 {
                #expect(abs(w - omegaMax) / omegaMax < 1e-12, "BZ 边界达带宽")
            }
        }
        // 零模：m = 0 → ω = 0（平移不变性/声学支 Goldstone 模）
        let iZero = ks.firstIndex { $0 == 0 } ?? -1
        #expect(iZero >= 0 && omegas[iZero] == 0)
        // 对称色散 ω(−k) = ω(k)
        for i in ks.indices { #expect(abs(omegas[i] - QFTBasicsMath.dispersion(-ks[i], a: a, v: v)) < 1e-18) }
    }

    @Test("物理律：零点能密度 u ∝ 1/a²（紫外截断敏感）——u·a² 截断扫描下严格守恒")
    func zeroPointScalingLaw() throws {
        let c = try constants()
        let hbar = try c.value("hbar")
        let v = 1e3
        let aList = [5e-10, 5e-11, 5e-12, 5e-13]
        var ua2Values: [Double] = []
        for a in aList {
            let e0 = QFTBasicsMath.zeroPointEnergy(N: 50, a: a, v: v, hbar: hbar)
            #expect(e0 > 0, "零点能恒正（真空不空）")
            let u = e0 / (50.0 * a)
            ua2Values.append(u * a * a)
        }
        for value in ua2Values.dropFirst() {
            #expect(abs(value / ua2Values[0] - 1) < 1e-12,
                    "u·a² = 常数：a → a/10 ⇒ u × 100（1/a² 发散的来源）")
        }
        // E₀ = Σℏω/2 与逐模求和一致
        let a = 5e-10
        let ks = QFTBasicsMath.waveVectors(N: 50, a: a)
        let direct = ks.reduce(0.0) { $0 + 0.5 * hbar * QFTBasicsMath.dispersion($1, a: a, v: v) }
        #expect(abs(direct - QFTBasicsMath.zeroPointEnergy(N: 50, a: a, v: v, hbar: hbar)) < 1e-30)
    }

    @Test("物理律：BE 占据 ⟨n_k⟩ 随 |k| 单调降、随 T 单调升；高温经典极限 ⟨n⟩ ≈ kT/ℏω − 1/2")
    func boseOccupancyLaw() throws {
        let c = try constants()
        let hbar = try c.value("hbar")
        let kB = try c.value("kB")
        let a = 5e-10, v = 1e3
        let ks = QFTBasicsMath.waveVectors(N: 50, a: a)
        let omegas = ks.map { QFTBasicsMath.dispersion($0, a: a, v: v) }

        // 零模占据 = 0（无热激发；脚本同口径特判）
        #expect(QFTBasicsMath.boseOccupancy(0, temperature: 300, hbar: hbar, kB: kB) == 0)

        // |k| 单调降（同一 T 下长波模优先占据；零模已特判 0，排除后检查）
        let occ = zip(ks, omegas).filter { $1 > 0 }.map { (k, w) in
            (abs(k), QFTBasicsMath.boseOccupancy(w, temperature: 10, hbar: hbar, kB: kB))
        }.sorted { $0.0 < $1.0 }
        for i in occ.indices.dropFirst() {
            #expect(occ[i].1 <= occ[i - 1].1 + 1e-12, "⟨n⟩ 随 |k| 递减")
        }
        // 温度单调升
        for w in omegas where w > 0 {
            let n10 = QFTBasicsMath.boseOccupancy(w, temperature: 10, hbar: hbar, kB: kB)
            let n300 = QFTBasicsMath.boseOccupancy(w, temperature: 300, hbar: hbar, kB: kB)
            #expect(n300 > n10)
            #expect(n10 >= 0 && n300 >= 0)
        }
        // 经典极限展开：x ≪ 1 时 1/(e^x−1) ≈ 1/x − 1/2 + x/12
        let wEdge = omegas.max() ?? 1   // BZ 边界模在 300 K 下 x ≈ 0.10
        let x = hbar * wEdge / (kB * 300)
        #expect(x < 0.15, "x = \(x) ≪ 1（经典区）")
        let exact = QFTBasicsMath.boseOccupancy(wEdge, temperature: 300, hbar: hbar, kB: kB)
        let classical = kB * 300 / (hbar * wEdge) - 0.5
        #expect(abs(exact - classical) / exact < 1e-3,
                "高温经典极限：\(exact) vs \(classical)")
    }

    // MARK: - fixtures 对拍

    @Test("fixtures：ideal exponent −4 参考线（主曲线无标签被 runner 过滤，幂律由物理律覆盖）")
    func casimirFixture() throws {
        let fx = try ModuleFixture.load("量子场论的基本图像_Casimir__v2022")
        let (x, y) = try fx.line(0, 1, label: "ideal exponent")
        #expect(y[0] == -4 && y[1] == -4, "理想指数参考线恒 −4")
        #expect(x[0] == 0 && x[1] == 1)
    }

    @Test("fixtures：BE 占据 T = 10 K / 300 K 双线各 50 点（< 1e-8，按数组序对拍）")
    func boseOccupancyMatchesFixtures() throws {
        let c = try constants()
        let hbar = try c.value("hbar")
        let kB = try c.value("kB")
        let a = 5e-10, v = 1e3     // 脚本固定 N = 50、a = 0.5 nm、v = 1 km/s

        let fx = try ModuleFixture.load("量子场论的基本图像_模式量子化__v2022")
        let ks = QFTBasicsMath.waveVectors(N: 50, a: a)
        let omegas = ks.map { QFTBasicsMath.dispersion($0, a: a, v: v) }

        let (_, y10) = try fx.line(0, 2, label: "T = 10 K")
        let (_, y300) = try fx.line(0, 2, label: "T = 300 K")
        #expect(y10.count == 50 && y300.count == 50)

        expectPointwiseClose(
            omegas.map { QFTBasicsMath.boseOccupancy($0, temperature: 10, hbar: hbar, kB: kB) },
            y10, tolerance: 1e-8, "BE 占据 T = 10 K")
        expectPointwiseClose(
            omegas.map { QFTBasicsMath.boseOccupancy($0, temperature: 300, hbar: hbar, kB: kB) },
            y300, tolerance: 1e-8, "BE 占据 T = 300 K")
        // fixture x 轴为脚本口径 k/a（含 a 量纲），与 k·a/(2π)（模块无量纲轴）线性对应：
        // 抽查首点确认同一 m 序（m = −25 起）
        let (xF, _) = try fx.line(0, 2, label: "T = 10 K")
        #expect(abs(xF[0] - ks[0] / a) / abs(xF[0]) < 1e-8,
                "x[0] = k[0]/a = \(xF[0])（脚本横轴口径）")
    }

    // MARK: - compute 结构与预算

    @Test("compute 结构：Casimir 双图（200/199 点）+ 摘要 + 理论卡；实时档预算")
    func casimirComputeStructure() async throws {
        let module = CasimirModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 2)
        for (index, chart) in result.charts.enumerated() {
            guard case .lineSeries(let c) = chart else {
                Issue.record("charts[\(index)] 应为 lineSeries"); return
            }
            #expect(c.series[0].points.count == (index == 0 ? 200 : 199))
        }
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 4)
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }

    @Test("compute 结构：模式量子化三图（色散/零点能/占据）+ 摘要 + 理论卡；实时档预算")
    func modeQuantizationComputeStructure() async throws {
        let module = ModeQuantizationModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 3)
        for (index, chart) in result.charts.enumerated() {
            guard case .lineSeries(let c) = chart else {
                Issue.record("charts[\(index)] 应为 lineSeries"); return
            }
            #expect(c.series[0].points.count == 50, "N = 50 模式")
        }
        #expect(result.summary.count == 5)
        #expect(result.theory?.formulas.count == 4)
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 16)
    }
}
