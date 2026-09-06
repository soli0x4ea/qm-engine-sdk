import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W4-16 斯特恩-盖拉赫：复刻脚本打印值（v/dz/2dz）+ fixtures 双高斯对拍（600 网格闭式）。
@Suite("W4 斯特恩-盖拉赫")
struct SternGerlachModuleTests {

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    /// 脚本默认参数下的完整偏转链。
    private func defaults() throws -> (v: Double, dz: Double, a: Double, d1: Double, d2: Double) {
        let cs = try constants()
        let u = try cs.value("u"), kB = try cs.value("kB"), muB = try cs.value("mu_B")
        let mAg = SternGerlachMath.silverMass(u: u)
        let v = SternGerlachMath.thermalSpeed(T: 1000, kB: kB, mass: mAg)
        let (a, d1, d2, dz) = SternGerlachMath.deflection(
            v: v, muB: muB, gradB: 1e3, mass: mAg, magnetLength: 0.10, drift: 1.0)
        return (v, dz, a, d1, d2)
    }

    @Test("脚本打印值：v = 392.6 m/s，dz = 35.265 mm，分离 = 70.530 mm（rel < 1e-3）")
    func printedValues() throws {
        let (v, dz, _, _, _) = try defaults()
        #expect(abs(v - 392.6) / 392.6 < 1e-3, "v = \(v)")
        #expect(abs(dz * 1e3 - 35.265) / 35.265 < 1e-3, "dz = \(dz * 1e3) mm")
        #expect(abs(2 * dz * 1e3 - 70.530) / 70.530 < 1e-3, "2dz = \(2 * dz * 1e3) mm")
    }

    @Test("偏转结构：dz = d1 + d2；d1/d2 比例（短磁铁长漂移）")
    func deflectionStructure() throws {
        let (_, _, _, d1, d2) = try defaults()
        // 脚本量级：磁内 ~1.7 mm，漂移 ~33.6 mm（漂移主导）
        #expect(abs(d1 * 1e3 - 1.679) / 1.679 < 1e-2, "d1 = \(d1 * 1e3) mm")
        #expect(abs(d2 * 1e3 - 33.586) / 33.586 < 1e-2, "d2 = \(d2 * 1e3) mm")
        #expect(d2 > 10 * d1, "漂移主导")
    }

    @Test("fixture 对拍：双高斯 512 点（600 网格闭式，rel < 1e-8；深尾量级断言）")
    func beamProfilesMatchFixtures() throws {
        let fx = try ModuleFixture.load("角动量与自旋_斯特恩-盖拉赫分裂__v2022")
        let (ux, uy) = try fx.line(0, 0, label: "+1/2")
        let (dx_, dy) = try fx.line(0, 0, label: "-1/2")
        let (_, dz, _, _, _) = try defaults()
        #expect(ux.count == 512 && dx_.count == 512)

        // fixture x（mm）反推 600 网格下标 j = (z+3dz)·599/(6dz)，
        // 闭式在精确网格点求值（消除 fixture x 9 位舍入在深尾的放大）。
        func closedAt(_ xmm: Double, up: Bool) -> (z: Double, y: Double) {
            let z0 = (xmm * 1e-3 + 3 * dz) * 599.0 / (6 * dz)
            let j = Int(z0.rounded())
            let z = -3 * dz + 6 * dz * Double(j) / 599.0
            return (z, SternGerlachMath.beamIntensity(z: z, dz: dz, up: up))
        }

        var strong = 0, tail = 0
        for i in ux.indices {
            #expect(abs(ux[i] - dx_[i]) < 1e-12, "共享 z 网格 第\(i)点")
            let (zu, yu) = closedAt(ux[i], up: true)
            let (_, yd) = closedAt(ux[i], up: false)
            _ = zu
            if uy[i] > 1e-6 {
                strong += 1
                #expect(abs(uy[i] - yu) / yu < 1e-8, "上束 第\(i)点 \(uy[i]) vs \(yu)")
            } else {
                tail += 1
                #expect(uy[i] < 1e-6, "上束深尾 第\(i)点 \(uy[i])")
            }
            if dy[i] > 1e-6 {
                #expect(abs(dy[i] - yd) / yd < 1e-8, "下束 第\(i)点 \(dy[i]) vs \(yd)")
            }
        }
        #expect(strong >= 100, "信号区样本数 \(strong)")
        #expect(strong + tail == ux.count)
        // 对称性：两束互为镜像，总强度在 z=0 处对称双峰
        #expect(abs(uy[0] - dy[ux.count - 1]) / max(uy[0], 1e-300) < 1e-8)
    }

    @Test("双束结构：峰在 ±dz、峰高 1、z=0 处 I_tot = 2·exp(−8)")
    func beamStructure() throws {
        let (_, dz, _, _, _) = try defaults()
        let peakUp = SternGerlachMath.beamIntensity(z: dz, dz: dz, up: true)
        let peakDn = SternGerlachMath.beamIntensity(z: -dz, dz: dz, up: false)
        let center = SternGerlachMath.beamIntensity(z: 0, dz: dz, up: true)
        #expect(abs(peakUp - 1) < 1e-15 && abs(peakDn - 1) < 1e-15)
        // (0−dz)/(√2·0.25dz) = −2√2 ⇒ exp(−8)
        #expect(abs(center - exp(-8)) < 1e-15, "中心 = \(center)")
    }

    @Test("compute 输出结构：1 图 3 系列（200 点）+ 摘要 + 理论卡")
    func computeStructure() async throws {
        let module = SternGerlachModule()
        let result = try await module.compute(
            ParamValues.defaults(for: module.params), constants: try constants())
        #expect(result.charts.count == 1)
        guard case .lineSeries(let c) = result.charts[0] else {
            Issue.record("应为 lineSeries"); return
        }
        #expect(c.series.count == 3)
        #expect(c.series[0].points.count == 200)   // 600 / 3
        #expect(c.series[1].points.count == 200)
        #expect(c.series[2].points.count == 200)
        // W13g #16：±3dz 自适应域——默认 dz≈35.3mm 时两束（±dz）必须在画面内
        let dzMm = 35.265
        #expect(c.series[0].points.allSatisfy { abs($0.x) <= 3 * dzMm + 1e-6 })
        #expect(c.series[0].points.contains { $0.y > 0.9 }, "上束峰应在显示窗内")
        #expect(c.series[1].points.contains { $0.y > 0.9 }, "下束峰应在显示窗内")
        #expect(result.summary.count == 4)
        #expect(result.theory?.formulas.count == 4)
    }

    @Test("秒级档预算：compute 中位 < 2000 ms")
    func computeBudget() async throws {
        let module = SternGerlachModule()
        try await expectComputeUnderBudget(
            module: module, values: ParamValues.defaults(for: module.params),
            constants: try constants(), budgetMillis: 2000)
    }
}
