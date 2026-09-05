import Foundation

/// 模块所属的 10 大分类（方案 §3.1 信息架构，与 47 篇系统笔记的章节对应）。
/// rawValue 即目录页显示名；W2 契约冻结后仅追加不改。
public enum ModuleCategory: String, Sendable, CaseIterable, Codable {
    case oldQuantum = "旧量子论与波粒二象性"
    case quantumStates = "量子态与数学结构"
    case schrodinger1D = "薛定谔方程与一维问题"
    case angularCentral = "角动量与中心力场"
    case formalTheory = "形式理论"
    case manyBody = "多体与凝聚态"
    case quantumInfo = "量子信息"
    case quantumOptics = "量子光学"
    case relativisticQFT = "相对论量子与场论"
    case gemology = "宝石学量子专题"
}

/// 交互性能档（调度策略与徽章颜色的唯一依据）。
/// - realtime: 拖动即重算，预算 16 ms（compute + render）
/// - seconds: 松手触发后台计算，预算 2 s，带轻量进度
/// - strategy: 策略档串行队列，首算 10 s、缓存命中 100 ms
public enum ComputeTier: String, Sendable, CaseIterable, Codable {
    case realtime
    case seconds
    case strategy

    /// 目录/模块页徽章文案
    public var badge: String {
        switch self {
        case .realtime: return "实时档"
        case .seconds: return "秒级档"
        case .strategy: return "策略档"
        }
    }
}

/// 学习难度二分（目录页筛选 chip：基础 / 专题）。
public enum Difficulty: String, Sendable, CaseIterable, Codable {
    case basic       // 基础
    case advanced    // 专题

    public var badge: String { self == .basic ? "基础" : "专题" }
}

/// 模块元数据——目录页、搜索、徽章、路由的唯一静态来源。
/// 63 个模块各持有一份；字段与开发计划 §3.2 契约一一对应。
public struct ModuleMeta: Sendable, Equatable {
    /// 全局唯一 id（建议格式：笔记编号_短名，如 "03_blackbody"）
    public let id: String
    /// 模块标题（如"黑体辐射 · 三律对比"）
    public let title: String
    /// 一句话简介（目录页副行）
    public let subtitle: String
    public let category: ModuleCategory
    /// 所属系统笔记编号（01–47）
    public let noteNumber: Int
    public let tier: ComputeTier
    public let difficulty: Difficulty
    /// 搜索关键词（物理量名、别名；标题与分类自动纳入搜索域）

    public let keywords: [String]

    public init(id: String, title: String, subtitle: String,
                category: ModuleCategory, noteNumber: Int,
                tier: ComputeTier, difficulty: Difficulty,
                keywords: [String] = []) {
        self.id = id
        self.title = title
        self.subtitle = subtitle
        self.category = category
        self.noteNumber = noteNumber
        self.tier = tier
        self.difficulty = difficulty
        self.keywords = keywords
    }

    /// 搜索命中：标题 / 副题 / 分类名 / 关键词 / 笔记编号（"03" 或 "3"）
    public func matches(query: String) -> Bool {
        let q = query.trimmingCharacters(in: .whitespaces)
        guard !q.isEmpty else { return true }
        let haystacks: [String] = [title, subtitle, category.rawValue] + keywords
        let lower = q.lowercased()
        if haystacks.contains(where: { $0.lowercased().contains(lower) }) { return true }
        if String(format: "%02d", noteNumber).contains(q) { return true }
        if let n = Int(q), n == noteNumber { return true }
        return false
    }
}
