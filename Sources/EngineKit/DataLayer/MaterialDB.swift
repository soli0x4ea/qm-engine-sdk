import Foundation

// MARK: - 材料数据行

/// 晶体场致色矿物（笔记 43《过渡金属离子致色的量子机制》）。
/// delta_cm 为八面体晶体场分裂 Δ_o（cm⁻¹），来源 Burns 1993 / Nassau 1983。
public struct CrystalFieldMineral: Sendable, Codable, Equatable, Identifiable {
    public let id: String
    public let name: String
    public let short: String
    /// 致色离子（如 "Cr3+"）
    public let ion: String
    /// 宿主晶格（如 "Al2O3"）
    public let host: String
    /// 八面体晶体场分裂（cm⁻¹）
    public let deltaCm: Double

    enum CodingKeys: String, CodingKey {
        case id, name, short, ion, host
        case deltaCm = "delta_cm"
    }
}

/// F 心卤化物（笔记 44《色心与晶格缺陷的量子描述》）。
/// 数据为文献汇编值（Fowler 1968；Klick & Schulman 1957）。
public struct FCenterHalide: Sendable, Codable, Equatable, Identifiable {
    public let id: String
    public let name: String
    /// 高频（光学）介电常数 ε∞ ≈ n_opt²
    public let epsInf: Double
    /// 实验 F 吸收带峰值能量（eV）
    public let fBandExpEV: Double
    /// 晶格常数（Å），仅用于 Mollwo-Ivey 经验律对照
    public let latticeA: Double

    enum CodingKeys: String, CodingKey {
        case id, name
        case epsInf = "eps_inf"
        case fBandExpEV = "f_band_exp_eV"
        case latticeA = "lattice_A"
    }
}

/// 拉曼振动模式（笔记 45《光谱学仪器的量子基础》）。
/// 1332 cm⁻¹ 为金刚石一阶拉曼位移（宝石学鉴别钻石与仿制品的公认实验值）。
public struct RamanMode: Sendable, Codable, Equatable, Identifiable {
    public let id: String
    public let name: String
    /// 振动波数（cm⁻¹）
    public let nuVcm: Double
    /// 洛伦兹半高全宽（cm⁻¹）
    public let fwhmCm: Double

    enum CodingKeys: String, CodingKey {
        case id, name
        case nuVcm = "nu_v_cm"
        case fwhmCm = "fwhm_cm"
    }
}

// MARK: - MaterialDB v1

/// MaterialDB v1（W7 基建）：晶体场 / F 心 / 拉曼 三张材料参数表。
///
/// JSON 随 EngineKit 打包（`Resources/materials_v1.json`），由
/// `tools/export_materials.py` 从 code/ 下 Python 源脚本静态解析导出；
/// 单测逐行比对随包资源与 `fixtures/materials/materials_v1.json`，
/// 保证数据层不被手改漂移（计划 §W7 基建条）。
public struct MaterialDB: Sendable, Equatable, Decodable {
    public static let schemaVersion = "materialdb/v1"

    public let schema: String
    public let generated: String
    /// 源脚本与指纹（16 位 sha256 前缀）
    public let provenance: [String: MaterialProvenance]
    public let crystalFieldMinerals: [CrystalFieldMineral]
    public let fCenterHalides: [FCenterHalide]
    public let ramanModes: [RamanMode]

    public struct MaterialProvenance: Sendable, Codable, Equatable {
        public let script: String
        public let sha256: String

        enum CodingKeys: String, CodingKey {
            case script
            case sha256 = "sha256_16"
        }
    }

    enum CodingKeys: String, CodingKey {
        case schema, generated, provenance
        case crystalFieldMinerals, fCenterHalides, ramanModes
    }

    public init(schema: String, generated: String,
                provenance: [String: MaterialProvenance],
                crystalFieldMinerals: [CrystalFieldMineral],
                fCenterHalides: [FCenterHalide],
                ramanModes: [RamanMode]) {
        self.schema = schema
        self.generated = generated
        self.provenance = provenance
        self.crystalFieldMinerals = crystalFieldMinerals
        self.fCenterHalides = fCenterHalides
        self.ramanModes = ramanModes
    }

    // MARK: 查询

    public func mineral(id: String) -> CrystalFieldMineral? {
        crystalFieldMinerals.first { $0.id == id }
    }

    public func halide(id: String) -> FCenterHalide? {
        fCenterHalides.first { $0.id == id }
    }

    public func ramanMode(id: String) -> RamanMode? {
        ramanModes.first { $0.id == id }
    }
}

public enum MaterialDBError: Error, Equatable, Sendable, LocalizedError {
    case resourceMissing
    case decodeFailed(String)
    case schemaMismatch(expected: String, actual: String)

    public var errorDescription: String? {
        switch self {
        case .resourceMissing:
            return "MaterialDB 随包资源 materials_v1.json 缺失"
        case .decodeFailed(let why):
            return "MaterialDB 解码失败：\(why)"
        case .schemaMismatch(let expected, let actual):
            return "MaterialDB 模式不匹配：期望 \(expected)，实际 \(actual)"
        }
    }
}

extension MaterialDB {
    /// 从随包资源加载（同名资源在 App 与测试中均为 Bundle.module）。
    /// 进程内缓存——MaterialDB 为不可变值类型，缓存安全。
    public static func load() throws -> MaterialDB { try cachedDB.load() }

    /// 从任意 URL 加载（fixtures 对拍 / 测试注入用）。
    public static func load(url: URL) throws -> MaterialDB {
        let data = try Data(contentsOf: url)
        let db: MaterialDB
        do {
            db = try JSONDecoder().decode(MaterialDB.self, from: data)
        } catch {
            throw MaterialDBError.decodeFailed("\(error)")
        }
        guard db.schema == schemaVersion else {
            throw MaterialDBError.schemaMismatch(expected: schemaVersion, actual: db.schema)
        }
        return db
    }
}

private final class MaterialDBCache: @unchecked Sendable {
    private let lock = NSLock()
    private var store: MaterialDB?

    func load() throws -> MaterialDB {
        lock.lock()
        defer { lock.unlock() }
        if let hit = store { return hit }
        guard let url = Bundle.module.url(forResource: "materials_v1", withExtension: "json")
        else { throw MaterialDBError.resourceMissing }
        let db = try MaterialDB.load(url: url)
        store = db
        return db
    }
}

private let cachedDB = MaterialDBCache()
