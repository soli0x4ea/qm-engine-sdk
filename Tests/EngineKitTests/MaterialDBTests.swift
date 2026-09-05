import Testing
import Foundation
@testable import EngineKit

/// MaterialDB v1（W7 基建）：随包资源 ↔ fixtures 逐行对拍 + 数据完整性。
/// fixtures 由 tools/export_materials.py 从 code/ 下 Python 源脚本静态解析导出，
/// 保证 Swift 侧材料参数表不被手改漂移（计划 §W7 基建条）。
@Suite("MaterialDB v1")
struct MaterialDBTests {

    private func bundled() throws -> MaterialDB { try MaterialDB.load() }

    /// fixtures 副本（Tests/EngineKitTests/Fixtures/materials_v1.json，随包 .copy）
    private func fixture() throws -> MaterialDB {
        let url = try #require(
            Bundle.module.url(forResource: "materials_v1", withExtension: "json",
                              subdirectory: "Fixtures"))
        return try MaterialDB.load(url: url)
    }

    // MARK: - 装载与模式

    @Test("随包资源可加载且模式为 materialdb/v1")
    func schemaAndLoading() throws {
        let db = try bundled()
        #expect(db.schema == MaterialDB.schemaVersion)
        #expect(!db.generated.isEmpty)
        // 三张表的源脚本指纹均登记（可追溯到 Python 源脚本）
        #expect(db.provenance["crystalField"]?.script == "43_过渡金属离子致色的量子机制_模型.py")
        #expect(db.provenance["fCenter"]?.script == "44_色心与晶格缺陷的量子描述_F心类氢模型.py")
        #expect(db.provenance["raman"]?.script == "45_光谱学仪器的量子基础_拉曼谱.py")
        for (key, prov) in db.provenance {
            #expect(prov.sha256.count == 16, "\(key) 源脚本指纹应为 16 位十六进制")
        }
    }

    @Test("三表规模：矿物 3 · 卤化物 6 · 拉曼模式 3")
    func tableSizes() throws {
        let db = try bundled()
        #expect(db.crystalFieldMinerals.count == 3)
        #expect(db.fCenterHalides.count == 6)
        #expect(db.ramanModes.count == 3)
    }

    // MARK: - 逐行 round-trip（模块用到的每一行）

    @Test("数据完整性：三表逐行与 fixtures 完全一致（逐位）")
    func everyRowRoundTrips() throws {
        let db = try bundled()
        let fx = try fixture()
        #expect(db.schema == fx.schema)
        #expect(db.provenance == fx.provenance)
        #expect(db.crystalFieldMinerals == fx.crystalFieldMinerals,
                "晶体场矿物表与 fixtures 不一致")
        #expect(db.fCenterHalides == fx.fCenterHalides, "F 心卤化物表与 fixtures 不一致")
        #expect(db.ramanModes == fx.ramanModes, "拉曼模式表与 fixtures 不一致")
    }

    @Test("晶体场矿物：Ruby/Emerald/Peridot 的 Δ_o 与 λ=10⁷/Δ 逐位一致")
    func crystalFieldRows() throws {
        let db = try bundled()
        let ruby = try #require(db.mineral(id: "ruby"))
        let emerald = try #require(db.mineral(id: "emerald"))
        let peridot = try #require(db.mineral(id: "peridot"))
        #expect(ruby.deltaCm == 18000.0)
        #expect(emerald.deltaCm == 16500.0)
        #expect(peridot.deltaCm == 9524.0)
        #expect(ruby.ion == "Cr3+" && ruby.host == "Al2O3")
        #expect(emerald.ion == "Cr3+" && emerald.host == "Be3Al2Si6O18")
        #expect(peridot.ion == "Fe2+" && peridot.host == "olivine")
        // 反比律在表内成立：λ·Δ = 10⁷ nm·cm⁻¹（笔记 43 式）
        for m in db.crystalFieldMinerals {
            #expect(abs((1.0e7 / m.deltaCm) * m.deltaCm - 1.0e7) < 1e-6)
        }
        // Δ_o 递减序：Ruby > Emerald > Peridot（对应 λ 递增，颜色由红到绿）
        #expect(ruby.deltaCm > emerald.deltaCm && emerald.deltaCm > peridot.deltaCm)
    }

    @Test("F 心卤化物：六种盐的 ε∞ / E_F(exp) / 晶格常数逐位一致，且 a 升序")
    func fCenterRows() throws {
        let db = try bundled()
        let expected: [(String, Double, Double, Double)] = [
            ("LiF", 1.93, 5.08, 4.03),
            ("NaF", 1.74, 3.70, 4.63),
            ("NaCl", 2.38, 2.75, 5.64),
            ("KCl", 2.22, 2.30, 6.29),
            ("KBr", 2.43, 2.06, 6.60),
            ("KI", 2.65, 1.87, 7.07),
        ]
        #expect(db.fCenterHalides.map(\.id) == expected.map(\.0))
        for (i, row) in db.fCenterHalides.enumerated() {
            let (name, eps, eExp, a) = expected[i]
            #expect(row.name == name)
            #expect(row.epsInf == eps, "\(name) ε∞")
            #expect(row.fBandExpEV == eExp, "\(name) E_F(exp)")
            #expect(row.latticeA == a, "\(name) 晶格常数")
        }
        // 晶格常数严格升序（Mollwo-Ivey 律的自变量）
        for i in 1..<db.fCenterHalides.count {
            #expect(db.fCenterHalides[i].latticeA > db.fCenterHalides[i - 1].latticeA)
        }
    }

    @Test("拉曼模式：金刚石 1332 cm⁻¹ 在表内，宽度全正")
    func ramanRows() throws {
        let db = try bundled()
        let diamond = try #require(db.ramanMode(id: "diamond_raman"))
        #expect(diamond.nuVcm == 1332.0, "金刚石一阶拉曼位移（宝石学公认值）")
        #expect(diamond.fwhmCm == 6.0)
        #expect(db.ramanModes.map(\.nuVcm) == [520.0, 1332.0, 2900.0])
        #expect(db.ramanModes.map(\.fwhmCm) == [10.0, 6.0, 14.0])
        for m in db.ramanModes {
            #expect(m.nuVcm > 0 && m.fwhmCm > 0, "\(m.id) 波数与线宽须为正")
        }
    }

    // MARK: - 错误处理

    @Test("模式不匹配即抛错，不静默降级")
    func schemaValidation() throws {
        // 构造一个模式串不同的 JSON，确认解码器拒绝
        let bogus = """
        {"schema":"materialdb/v9","generated":"","provenance":{},
         "crystalFieldMinerals":[],"fCenterHalides":[],"ramanModes":[]}
        """
        let url = FileManager.default.temporaryDirectory
            .appendingPathComponent("materialdb_bogus_\(UUID().uuidString).json")
        try bogus.write(to: url, atomically: true, encoding: .utf8)
        defer { try? FileManager.default.removeItem(at: url) }
        #expect(throws: MaterialDBError.self) { try MaterialDB.load(url: url) }
    }
}
