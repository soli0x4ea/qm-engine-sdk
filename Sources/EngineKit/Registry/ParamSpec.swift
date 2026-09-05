import Foundation

// MARK: - 刻度

/// 滑杆与坐标轴的刻度类型。
public enum ParamScale: String, Sendable, Codable {
    case linear
    case log
}

// MARK: - 四类控件规格（开发计划 §3.2 契约）

/// ① 连续滑杆（线性 / 对数）。
/// 值域与默认值均为**显示单位**；computeUnit 非空时，引擎在调用 compute 前
/// 经 UnitKit 换算到计算单位（方案 §5.4：计算按模块声明单位制，换算只在边界）。
public struct SliderSpec: Sendable, Equatable {
    public let key: String
    public let title: String
    /// 物理量符号（如 "T"、"λ"，用于标签排版）
    public let symbol: String?
    /// 显示单位（如 "K"、"nm"）
    public let unit: String
    /// 显示单位下的取值范围（对数刻度要求全正）
    public let range: ClosedRange<Double>
    public let defaultValue: Double
    public let scale: ParamScale
    /// 步进（nil = 连续）
    public let step: Double?
    /// 显示保留小数位
    public let decimalPlaces: Int
    /// 计算单位（nil = 与显示单位一致）
    public let computeUnit: String?

    public init(key: String, title: String, symbol: String? = nil, unit: String,
                range: ClosedRange<Double>, defaultValue: Double,
                scale: ParamScale = .linear, step: Double? = nil,
                decimalPlaces: Int = 2, computeUnit: String? = nil) {
        precondition(range.contains(defaultValue), "SliderSpec[\(key)]: 默认值须在范围内")
        if scale == .log {
            precondition(range.lowerBound > 0 && defaultValue > 0, "SliderSpec[\(key)]: 对数刻度要求正值")
        }
        self.key = key
        self.title = title
        self.symbol = symbol
        self.unit = unit
        self.range = range
        self.defaultValue = defaultValue
        self.scale = scale
        self.step = step
        self.decimalPlaces = decimalPlaces
        self.computeUnit = computeUnit
    }
}

/// ② 离散分段选择的单个选项。
public struct DiscreteOption: Sendable, Equatable, Identifiable {
    public let id: String
    public let title: String
    /// 选项副行（如金属的功函数 "W = 4.5 eV"）
    public let subtitle: String?

    public init(id: String, title: String, subtitle: String? = nil) {
        self.id = id
        self.title = title
        self.subtitle = subtitle
    }
}

/// ② 离散分段控件（如光电效应的金属选择）。
public struct DiscreteSpec: Sendable, Equatable {
    public let key: String
    public let title: String
    public let options: [DiscreteOption]
    public let defaultOptionID: String

    public init(key: String, title: String, options: [DiscreteOption],
                defaultOptionID: String) {
        precondition(options.contains(where: { $0.id == defaultOptionID }),
                     "DiscreteSpec[\(key)]: defaultOptionID 不在选项中")
        self.key = key
        self.title = title
        self.options = options
        self.defaultOptionID = defaultOptionID
    }
}

/// ③ 多值对比组的候选项（如黑体辐射的"太阳 5772 K / 白炽灯 2800 K / CMB 2.725 K"）。
public struct CompareCandidate: Sendable, Equatable, Identifiable {
    public let id: String
    public let label: String
    /// 候选值（显示单位，语义由模块定义）
    public let value: Double

    public init(id: String, label: String, value: Double) {
        self.id = id
        self.label = label
        self.value = value
    }
}

/// ③ 多值对比组：多选若干候选值同图对比（曲线族叠加）。
public struct MultiCompareSpec: Sendable, Equatable {
    public let key: String
    public let title: String
    public let candidates: [CompareCandidate]
    public let defaultSelectionIDs: [String]
    /// 同时选中上限（图例可读性，默认 4）
    public let maxSelection: Int

    public init(key: String, title: String, candidates: [CompareCandidate],
                defaultSelectionIDs: [String], maxSelection: Int = 4) {
        precondition(!defaultSelectionIDs.isEmpty, "MultiCompareSpec[\(key)]: 至少默认选中一项")
        precondition(maxSelection >= 1 && defaultSelectionIDs.count <= maxSelection,
                     "MultiCompareSpec[\(key)]: 默认选中数超过上限")
        self.key = key
        self.title = title
        self.candidates = candidates
        self.defaultSelectionIDs = defaultSelectionIDs
        self.maxSelection = maxSelection
    }
}

/// ④ 只读常量卡：展示当前 ConstantsSet 中的 CODATA 值。
/// 值不参与 ParamValues 传递——compute 直接从 constants 参数取。
public struct ConstantSpec: Sendable, Equatable {
    public let key: String
    public let title: String
    /// 附加说明（如 "普朗克常量"）
    public let note: String?

    public init(key: String, title: String, note: String? = nil) {
        self.key = key
        self.title = title
        self.note = note
    }
}

/// 参数规格的四类封装——参数面板据此路由到对应控件。
public enum ParamSpec: Sendable, Equatable {
    case slider(SliderSpec)
    case discrete(DiscreteSpec)
    case multiCompare(MultiCompareSpec)
    case constant(ConstantSpec)

    public var key: String {
        switch self {
        case .slider(let s): return s.key
        case .discrete(let s): return s.key
        case .multiCompare(let s): return s.key
        case .constant(let s): return s.key
        }
    }
}

// MARK: - 参数值（用户输入的运行时集合，均为显示单位）

/// compute 的唯一输入。纯值类型，UI 每次交互生成新值。
public struct ParamValues: Sendable, Equatable {
    /// 滑杆值（显示单位）
    public var sliders: [String: Double]
    /// 离散选择（key → optionID）
    public var discretes: [String: String]
    /// 多值对比（key → 选中的 candidateID 集合）
    public var multiCompares: [String: Set<String>]

    public init(sliders: [String: Double] = [:],
                discretes: [String: String] = [:],
                multiCompares: [String: Set<String>] = [:]) {
        self.sliders = sliders
        self.discretes = discretes
        self.multiCompares = multiCompares
    }

    /// 从模块参数规格生成默认值（首次进入模块页 / 重置按钮共用）。
    public static func defaults(for params: [ParamSpec]) -> ParamValues {
        var v = ParamValues()
        for spec in params {
            switch spec {
            case .slider(let s):
                v.sliders[s.key] = s.defaultValue
            case .discrete(let s):
                v.discretes[s.key] = s.defaultOptionID
            case .multiCompare(let s):
                v.multiCompares[s.key] = Set(s.defaultSelectionIDs)
            case .constant:
                break  // 常量卡无输入态
            }
        }
        return v
    }

    // MARK: 便捷取值（缺省回退规格默认值由调用方保证；此处缺 key 直接致命错误，
    // 因为 UI 层永远以 defaults 为基底增量修改，正常路径不会缺）

    public func slider(_ key: String) -> Double {
        guard let v = sliders[key] else {
            fatalError("ParamValues: 缺少滑杆 \(key)（应以 defaults(for:) 为基底）")
        }
        return v
    }

    public func discrete(_ key: String) -> String {
        guard let v = discretes[key] else {
            fatalError("ParamValues: 缺少离散选择 \(key)")
        }
        return v
    }

    public func multiCompare(_ key: String) -> Set<String> {
        multiCompares[key] ?? []
    }

    /// 滑杆值换算到计算单位（引擎在 compute 入口统一调用）。
    /// computeUnit 为 nil 时原值返回。
    public func sliderComputeValue(_ spec: SliderSpec,
                                   constants: ConstantsSet) throws -> Double {
        let display = slider(spec.key)
        guard let target = spec.computeUnit else { return display }
        return try UnitKit.convert(display, from: spec.unit, to: target,
                                   constants: constants)
    }
}
