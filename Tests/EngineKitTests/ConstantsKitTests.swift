import Testing
import Foundation
@testable import EngineKit

/// CODATA 常量库验收：与 Python 端 constants.py 的 C 字典逐位一致。
/// 参考值由 tools/gen_reference_fixtures.py 从 data/constants_*.json 导出
/// （Tests/EngineKitTests/Fixtures/constants_reference.json）。
@Suite("ConstantsKit")
struct ConstantsKitTests {

    @Test("双版本加载，各 35 键")
    func loadBothVersions() throws {
        for v in CODATAVersion.allCases {
            let set = try ConstantsSet.load(v)
            #expect(set.count == 35, "\(v) 应含 35 个常量")
            #expect(set.meta.codataVersion == v.rawValue)
        }
    }

    @Test("SI 定义精确值（2019 重定义后定义常量）")
    func exactConstants() throws {
        let s = try ConstantsSet.load(.v2018)
        #expect(try s.value("h") == 6.62607015e-34)
        #expect(try s.value("hbar") == 1.054571817e-34)
        #expect(try s.value("e") == 1.602176634e-19)
        #expect(try s.value("c") == 299_792_458.0)
        #expect(try s.value("kB") == 1.380649e-23)
        #expect(try s.value("NA") == 6.02214076e23)
    }

    @Test("测量常量携带不确定度")
    func uncertainConstants() throws {
        let s = try ConstantsSet.load(.v2018)
        let me = try #require(s["m_e"])
        #expect(me.value == 9.1093837015e-31)
        #expect(me.unc == 2.8e-40)
        #expect(me.exact == false)
        #expect(me.unit == "kg")
    }

    @Test("双版本差异：依赖测量的常量在 2022 调整中变化")
    func versionDrift() throws {
        let s18 = try ConstantsSet.load(.v2018)
        let s22 = try ConstantsSet.load(.v2022)
        let h18 = try s18.value("h"), h22 = try s22.value("h")
        let e18 = try s18.value("e"), e22 = try s22.value("e")
        let a18 = try s18.value("alpha"), a22 = try s22.value("alpha")
        let me18 = try s18.value("m_e"), me22 = try s22.value("m_e")
        #expect(h18 == h22)
        #expect(e18 == e22)
        #expect(a18 != a22)
        #expect(me18 != me22)
    }

    @Test("缺失键抛错，不静默回退")
    func missingKeyThrows() throws {
        let s = try ConstantsSet.load(.v2018)
        #expect(throws: ConstantsError.self) { try s.value("nonexistent_key") }
    }

    @Test("全量 35 键 × 2 版本与 Python 导出参考逐位一致")
    func fullTableAgainstPythonReference() throws {
        let ref = try Self.loadReference()
        for v in CODATAVersion.allCases {
            let set = try ConstantsSet.load(v)
            let expected = try #require(ref[v.rawValue] as? [String: Double])
            #expect(expected.count == set.count, "\(v) 键数一致")
            for (k, want) in expected {
                #expect(try set.value(k) == want, "\(v)/\(k)")
            }
        }
    }

    // MARK: - 参考值加载

    private static func loadReference() throws -> [String: Any] {
        let url = try #require(
            Bundle.module.url(forResource: "constants_reference", withExtension: "json",
                              subdirectory: "Fixtures")
        )
        let data = try Data(contentsOf: url)
        return try #require(try JSONSerialization.jsonObject(with: data) as? [String: Any])
    }
}
