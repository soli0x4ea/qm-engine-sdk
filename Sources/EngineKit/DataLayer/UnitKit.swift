import Foundation

/// 物理量纲（换算只在同量纲内进行）
public enum QuantityKind: String, Sendable, CaseIterable {
    case energy        // J 系
    case length        // m 系
    case time          // s 系
    case frequency     // Hz 系
    case waveNumber    // m^-1 / cm^-1
    case mass          // kg / u
}

public enum UnitKitError: Error, Equatable, Sendable {
    /// 未知单位符号
    case unknownUnit(String)
    /// 量纲不匹配（如 J → nm）
    case kindMismatch(from: String, to: String)
    /// 常量缺失（依赖 ConstantsSet 的换算，如 a0 / kB）
    case constantMissing(String)
}

/// 单位换算层（方案 5.4）：「计算单位 → 显示单位」单向管线。
///
/// W1 范围：常用物理量的静态换算表 + 依赖 CODATA 的换算（K↔J、a0、u、cm^-1）。
/// W2 将在此之上声明 ParamSpec 的刻度（线性 / 对数）与单位选择。
///
/// 设计约束：换算只发生在渲染边界，引擎计算一律使用模块声明的计算单位，
/// 保证与 Python 原版逐位一致。
public enum UnitKit {

    /// 单位符号（与模块元数据 / 显示层共用词表）
    public static let symbols: [QuantityKind: [String]] = [
        .energy: ["J", "eV", "meV", "keV", "MeV", "K", "Wh"],
        .length: ["m", "nm", "μm", "mm", "a0"],
        .time: ["s", "fs", "ps", "ns"],
        .frequency: ["Hz", "kHz", "MHz", "GHz", "THz"],
        .waveNumber: ["m^-1", "cm^-1"],
        .mass: ["kg", "u"],
    ]

    /// 换算系数：1 单位 = `siFactor` 个 SI 基本单位（J / m / s / Hz / m^-1 / kg）
    /// 不依赖常量的部分为纯静态表；依赖常量的（K、a0、u）在 `factor(_:to:constants:)` 中解析。
    static let staticFactors: [String: (kind: QuantityKind, si: Double)] = [
        "J": (.energy, 1),
        "eV": (.energy, 1.602176634e-19),           // exact（SI 定义）
        "meV": (.energy, 1.602176634e-22),
        "keV": (.energy, 1.602176634e-16),
        "MeV": (.energy, 1.602176634e-13),
        "m": (.length, 1),
        "nm": (.length, 1e-9),
        "μm": (.length, 1e-6),
        "mm": (.length, 1e-3),
        "s": (.time, 1),
        "fs": (.time, 1e-15),
        "ps": (.time, 1e-12),
        "ns": (.time, 1e-9),
        "Hz": (.frequency, 1),
        "kHz": (.frequency, 1e3),
        "MHz": (.frequency, 1e6),
        "GHz": (.frequency, 1e9),
        "THz": (.frequency, 1e12),
        "m^-1": (.waveNumber, 1),
        "cm^-1": (.waveNumber, 100),
        "kg": (.mass, 1),
    ]

    /// 依赖常量的单位 → (量纲, 取值键)
    static let constantUnits: [String: (kind: QuantityKind, key: String)] = [
        "K": (.energy, "kB"),        // 温度的能量当量：1 K ↔ kB J
        "a0": (.length, "a0"),       // 玻尔半径
        "u": (.mass, "u"),           // 原子质量单位
    ]

    /// 波数 → 能量需经 E = h c ν̃（能量与波数跨量纲，单独提供方法）
    public static func energyFromWaveNumber(_ nuTildePerMeter: Double,
                                            constants: ConstantsSet) throws -> Double {
        let h = try constants.value("h")
        let c = try constants.value("c")
        return h * c * nuTildePerMeter
    }

    /// 换算系数：`value_in_to = value_in_from × factor(from:from, to:to)`
    /// （from 单位的 SI 量 ÷ to 单位的 SI 量；如 factor(eV→J) = 1.602…e-19）
    public static func factor(from: String, to: String,
                              constants: ConstantsSet? = nil) throws -> Double {
        let f = try siFactor(from, constants: constants)
        let t = try siFactor(to, constants: constants)
        guard f.kind == t.kind else {
            throw UnitKitError.kindMismatch(from: from, to: to)
        }
        return f.si / t.si
    }

    /// 便捷换算：`convert(2.5, from: "eV", to: "meV") == 2500`
    public static func convert(_ value: Double, from: String, to: String,
                               constants: ConstantsSet? = nil) throws -> Double {
        try value * factor(from: from, to: to, constants: constants)
    }

    static func siFactor(_ unit: String, constants: ConstantsSet?) throws -> (kind: QuantityKind, si: Double) {
        if let s = staticFactors[unit] { return s }
        if let (kind, key) = constantUnits[unit] {
            guard let cs = constants else { throw UnitKitError.constantMissing(key) }
            return (kind, try cs.value(key))
        }
        throw UnitKitError.unknownUnit(unit)
    }
}
