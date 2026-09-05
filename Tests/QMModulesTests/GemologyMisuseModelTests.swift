import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// 笔记 47《宝石学中的量子误用辨析_模型》模块测试。
/// fixture 为双 3D Bloch 球 + 能级示意图（零线条数据）→ 结构级对拍；
/// 物理律（纯/混对照、R 线、MaterialDB round-trip）扛主力。
@Suite(.serialized)
struct GemologyMisuseModelTests {

    private let module = GemologyMisuseModelModule()

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    private func values(delta: Double = .pi / 2) -> ParamValues {
        var v = ParamValues.defaults(for: module.params)
        v.sliders["delta"] = delta
        return v
    }

    // MARK: fixture 结构对拍

    @Test("fixture：结构级对拍（双 3D Bloch 球 + 红宝石能级图，零线条）")
    func fixtureStructure() throws {
        let fixture = try ModuleFixture.load("47_宝石学中的量子误用辨析_模型__v2022")
        #expect(fixture.document["status"] as? String == "ok")
        let figures = try #require(fixture.document["figures"] as? [[String: Any]])
        #expect(figures.count == 2, "脚本产出图 1（双 Bloch 球）+ 图 2（红宝石能级）")
        let axes0 = try #require(figures[0]["axes"] as? [[String: Any]])
        #expect(axes0.count == 2)
        #expect((axes0[0]["title"] as? String)?.contains("PURE separable") == true)
        #expect((axes0[1]["title"] as? String)?.contains("maximally MIXED") == true)
        let axes1 = try #require(figures[1]["axes"] as? [[String: Any]])
        #expect((axes1[0]["title"] as? String)?.contains("Ruby Cr3+ level scheme") == true)
        for axes in [axes0, axes1] {
            for ax in axes {
                let lines = ax["lines"] as? [[String: Any]] ?? []
                #expect(lines.isEmpty, "3D 线框 / hlines 不导出 lines——对拍由物理律承担")
            }
        }
    }

    // MARK: 物理律 1：分离纯态 |r|=1、S=0、r=(cos δ, sin δ, 0)

    @Test("物理律1：双折射后单光子纯态 r=(cos δ, sin δ, 0)、|r|=1、S=0")
    func separablePureState() {
        for j in 0..<9 {
            let delta = Double(j) * .pi / 4
            let r = GemologyMisuseMath.separableBloch(delta: delta)
            #expect(abs(r.rx - cos(delta)) < 1e-12)
            #expect(abs(r.ry - sin(delta)) < 1e-12)
            #expect(abs(r.rz) < 1e-12, "等幅叠加 rz = p0−p1 = 0")
            let norm = sqrt(r.rx * r.rx + r.ry * r.ry + r.rz * r.rz)
            #expect(abs(norm - 1) < 1e-12, "球面纯态")
            #expect(abs(BlochMath.vonNeumannEntropyBits(r)) < 1e-12, "S = 0 bit")
        }
    }

    // MARK: 物理律 2：最大混态 r=0、S=1；径向差 = 1

    @Test("物理律2：约化最大混态 r=0、S=1 bit；纯/混径向差 = 1（双折射不增纠缠）")
    func mixedStateAndRadialContrast() {
        let rEnt = GemologyMisuseMath.maximallyMixedBloch()
        #expect(rEnt.rx == 0 && rEnt.ry == 0 && rEnt.rz == 0)
        #expect(abs(BlochMath.vonNeumannEntropyBits(rEnt) - 1) < 1e-12)
        #expect(abs(BlochMath.purity(rEnt) - 0.5) < 1e-12)
        for j in 0..<5 {
            let rSep = GemologyMisuseMath.separableBloch(delta: Double(j) * .pi / 2)
            let normSep = sqrt(rSep.rx * rSep.rx + rSep.ry * rSep.ry + rSep.rz * rSep.rz)
            #expect(abs(normSep - 1) < 1e-12)
            #expect(abs(normSep - 1) - 0 < 1e-12, "径向差 |r纯| − |r混| = 1 − 0 = 1")
        }
    }

    // MARK: 物理律 3：红宝石 R 线与 ZFS

    @Test("物理律3：R1 = 694.35 nm ↔ 脚本 694.3 nm；R1−R2 = 29 cm⁻¹；R1 能量 1.786 eV")
    func rubyRLines() {
        let r1Nm = GemologyMisuseMath.wavenumberToNm(GemologyMisuseMath.e2ER1)
        let r2Nm = GemologyMisuseMath.wavenumberToNm(GemologyMisuseMath.e2ER2)
        #expect(abs(r1Nm - 694.3) < 0.1, "脚本标注 694.3 nm")
        #expect(abs(r2Nm - 692.9) < 0.1, "脚本标注 692.9 nm")
        #expect(abs(GemologyMisuseMath.zfs2E - 29) < 1e-12, "R1/R2 零场分裂 ~29 cm⁻¹")
        let eR1 = GemologyMisuseMath.energyEV(wavelengthNm: 694.3)
        // 脚本 __main__ 回显 energy_ev(694.3)；hc = 1239.841984 eV·nm
        #expect(abs(eR1 - 1.785755) < 1e-4)
        #expect(abs(GemologyMisuseMath.energyEV(wavelengthNm: r1Nm) - eR1) < 2e-4,
                "波数路径与脚本波长路径一致（λ 差 0.047 nm → ΔE ≈ 1.2e-4 eV）")
        #expect(abs(1e7 / 555.6 - GemologyMisuseMath.e4T2) < 30,
                "4T2 ≈ 555 nm 吸收带（脚本注释口径）")
    }

    // MARK: 物理律 4：MaterialDB ruby round-trip + 三行核验

    @Test("物理律4：MaterialDB 三行逐行核验（ruby 18000 = E(4T2)；emerald/peridot 对源）")
    func materialDBRoundTrip() throws {
        let db = try MaterialDB.load()
        let ruby = try #require(db.mineral(id: "ruby"))
        #expect(ruby.ion == "Cr3+" && ruby.host == "Al2O3")
        #expect(ruby.deltaCm == 18000.0)
        #expect(ruby.deltaCm == GemologyMisuseMath.e4T2,
                "ruby Δ_o 与 4T2 宽吸收带同源（round-trip）")
        let emerald = try #require(db.mineral(id: "emerald"))
        #expect(emerald.ion == "Cr3+" && emerald.host == "Be3Al2Si6O18")
        #expect(emerald.deltaCm == 16500.0, "43 号脚本源值")
        let peridot = try #require(db.mineral(id: "peridot"))
        #expect(peridot.ion == "Fe2+" && peridot.host == "olivine")
        #expect(peridot.deltaCm == 9524.0, "43 号脚本源值（~1050 nm 近红外主带）")
        #expect(abs(1e7 / peridot.deltaCm - 1050) < 1)
    }

    // MARK: compute 输出与预算

    @Test("compute：能级图数据 + S(|r|) 曲线 + 实时档 16 ms 预算")
    func computeAndBudget() async throws {
        let cs = try constants()
        let result = try await module.compute(values(), constants: cs)
        #expect(result.charts.count == 2)
        // 能级图：5 能级升序、4 跃迁
        guard case .levelDiagram(let ld) = result.charts[1] else {
            Issue.record("第 2 图应为能级图")
            return
        }
        #expect(ld.levels.count == 5)
        #expect(ld.levels.map(\.energy) == ld.levels.map(\.energy).sorted(), "能量升序")
        #expect(ld.transitions.count == 4, "R1/R2 发射 + 4T2/4T1 吸收")
        #expect(ld.levels[0].energy == 0 && ld.levels[3].energy == 18000)
        // S(|r|) 曲线端点：S(0)=1、S(1)=0
        let chart = try #require(chartLineSeries(result, 0))
        let pts = chart.series[0].points
        #expect(abs(pts.first!.y - 1) < 1e-12)
        #expect(abs(pts.last!.y) < 1e-12)
        try await expectComputeUnderBudget(module: module, values: values(),
                                           constants: cs, budgetMillis: 16)
    }
}
