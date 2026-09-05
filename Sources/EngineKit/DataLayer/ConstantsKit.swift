import Foundation

/// CODATA 推荐值版本（2018 / 2022 双版本随包内置）
public enum CODATAVersion: String, Sendable, CaseIterable, Codable {
    case v2018 = "2018"
    case v2022 = "2022"
}

/// 单个物理常量：数值 + 单位 + 相对不确定度元数据
public struct PhysicalConstant: Sendable, Codable, Equatable {
    public let key: String
    public let name: String
    public let symbol: String
    public let value: Double
    public let unit: String
    public let exact: Bool
    public let unc: Double?
    public let note: String?

    enum CodingKeys: String, CodingKey {
        case name, symbol, value, unit, exact, unc, note
    }

    public init(key: String, name: String, symbol: String, value: Double,
                unit: String, exact: Bool, unc: Double? = nil, note: String? = nil) {
        self.key = key
        self.name = name
        self.symbol = symbol
        self.value = value
        self.unit = unit
        self.exact = exact
        self.unc = unc
        self.note = note
    }

    // Codable 只覆盖数据字段；key 由外层字典键提供，解码后由 ConstantsSet 补写
    public init(from decoder: Decoder) throws {
        let c = try decoder.container(keyedBy: CodingKeys.self)
        key = ""
        name = try c.decode(String.self, forKey: .name)
        symbol = try c.decode(String.self, forKey: .symbol)
        value = try c.decode(Double.self, forKey: .value)
        unit = try c.decode(String.self, forKey: .unit)
        // CODATA 源数据仅对 SI 定义常量标注 exact；测量/导出常量缺省为 false
        exact = try c.decodeIfPresent(Bool.self, forKey: .exact) ?? false
        unc = try c.decodeIfPresent(Double.self, forKey: .unc)
        note = try c.decodeIfPresent(String.self, forKey: .note)
    }

    public func encode(to encoder: Encoder) throws {
        var c = encoder.container(keyedBy: CodingKeys.self)
        try c.encode(name, forKey: .name)
        try c.encode(symbol, forKey: .symbol)
        try c.encode(value, forKey: .value)
        try c.encode(unit, forKey: .unit)
        try c.encode(exact, forKey: .exact)
        try c.encodeIfPresent(unc, forKey: .unc)
        try c.encodeIfPresent(note, forKey: .note)
    }
}

/// 常量集元数据（出处与核对说明）
public struct CODATAMeta: Sendable, Codable {
    public let codataVersion: String
    public let adjustment: String
    public let publication: String
    public let doi: String?
    public let sourceURL: String
    public let retrieved: String

    enum CodingKeys: String, CodingKey {
        case codataVersion = "codata_version"
        case adjustment, publication, doi
        case sourceURL = "source_url"
        case retrieved
    }
}

public enum ConstantsError: Error, Equatable, Sendable {
    /// 常量集中不存在该键（键名与 Python 端 constants.py 的 C 字典一致）
    case missingKey(String, CODATAVersion)
    /// 随包资源缺失或损坏
    case resourceBroken(CODATAVersion, String)
}

/// 一个 CODATA 版本的完整常量集。
///
/// 键名与 Python `constants.py` 的 `C` 字典完全一致（c / h / hbar / e / kB / m_e / a0 / …，共 35 键），
/// 保证移植模块 `C["m_e"]` ↔ `try constants.value("m_e")` 逐位一致。
public struct ConstantsSet: Sendable {
    public let version: CODATAVersion
    public let meta: CODATAMeta
    public let table: [String: PhysicalConstant]

    public subscript(key: String) -> PhysicalConstant? { table[key] }

    /// 数值访问（对齐 Python `C[key]`；缺失即抛错，不静默回退）
    public func value(_ key: String) throws -> Double {
        guard let c = table[key] else { throw ConstantsError.missingKey(key, version) }
        return c.value
    }

    /// 便捷访问：多个键一次取全（模块元数据声明依赖常量时使用）
    public func values(_ keys: [String]) throws -> [String: Double] {
        var out: [String: Double] = [:]
        for k in keys { out[k] = try value(k) }
        return out
    }

    public var count: Int { table.count }
}

private struct ConstantsFile: Codable {
    let meta: CODATAMeta
    let constants: [String: PhysicalConstant]
}

extension ConstantsSet {
    /// 从随包资源加载（W1 起双版本 JSON 随 EngineKit 打包）。
    /// 同一进程内版本级缓存——常量集为不可变值类型，缓存安全。
    public static func load(_ version: CODATAVersion) throws -> ConstantsSet {
        try cached.load(version)
    }
}

private final class ConstantsCache: @unchecked Sendable {
    private let lock = NSLock()
    private var store: [CODATAVersion: ConstantsSet] = [:]

    func load(_ version: CODATAVersion) throws -> ConstantsSet {
        lock.lock()
        defer { lock.unlock() }
        if let hit = store[version] { return hit }
        let set = try decode(version)
        store[version] = set
        return set
    }

    private func decode(_ version: CODATAVersion) throws -> ConstantsSet {
        let resource = "constants_\(version.rawValue)"
        guard let url = Bundle.module.url(forResource: resource, withExtension: "json"),
              let data = try? Data(contentsOf: url) else {
            throw ConstantsError.resourceBroken(version, "bundle resource \(resource).json not found")
        }
        do {
            let file = try JSONDecoder().decode(ConstantsFile.self, from: data)
            var table = file.constants
            for (k, v) in table { table[k] = PhysicalConstant(
                key: k, name: v.name, symbol: v.symbol, value: v.value,
                unit: v.unit, exact: v.exact, unc: v.unc, note: v.note) }
            return ConstantsSet(version: version, meta: file.meta, table: table)
        } catch {
            throw ConstantsError.resourceBroken(version, "decode failed: \(error)")
        }
    }
}

private let cached = ConstantsCache()
