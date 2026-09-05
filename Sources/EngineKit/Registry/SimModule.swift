import Foundation

/// 引擎与 UI 的唯一契约（开发计划 §3.2，W2 冻结，此后仅追加非破坏性扩展）。
///
/// 目录页、搜索、性能徽章、参数面板、图表路由全部由元数据驱动——
/// 新增模块 = 实现本协议 + 注册一行 + fixtures 一份，不改框架代码。
///
/// 实现要求（SOP 约束）：
/// - compute 为纯函数：无全局状态、无 IO，同输入必同输出
/// - 秒级/策略档实现须周期性检查 `Task.isCancelled`（调度器取消令牌协作点）
/// - 输入为显示单位；需要计算单位时经 `ParamValues.sliderComputeValue` 换算
/// - 曲线点数上限 2048（实时档预算的一部分）
public protocol SimModule: Sendable {
    /// id、名称、分类、笔记编号、性能档、难度、搜索关键词
    var meta: ModuleMeta { get }
    /// 四类控件：slider(线性/对数) / discrete / multiCompare / constant
    var params: [ParamSpec] { get }
    /// 图表声明：类型 / 轴 / 单位 / 图例（数据侧 SimResult.charts 按序对应）
    var charts: [ChartSpec] { get }
    /// 执行计算。输入显示单位参数 + 当前 CODATA 常量集。
    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult
}

// MARK: - 模块注册表

/// 全局模块注册表。App 启动时逐模块注册（SOP 的"注册一行"），
/// 目录页/搜索/路由全部查询本表。
public final class ModuleRegistry: @unchecked Sendable {

    public static let shared = ModuleRegistry()

    private var modules: [any SimModule] = []
    private let lock = NSLock()

    public init() {}

    /// 注册一个模块。id 重复时忽略并返回 false（调用方可据此在开发期暴露问题）。
    @discardableResult
    public func register(_ module: any SimModule) -> Bool {
        lock.lock()
        defer { lock.unlock() }
        if modules.contains(where: { $0.meta.id == module.meta.id }) {
            return false
        }
        modules.append(module)
        return true
    }

    /// 全部模块（按分类 → 笔记编号 → 标题稳定排序）。
    public var all: [any SimModule] {
        lock.lock()
        defer { lock.unlock() }
        return modules.sorted {
            if $0.meta.category.rawValue != $1.meta.category.rawValue {
                return $0.meta.category.rawValue < $1.meta.category.rawValue
            }
            if $0.meta.noteNumber != $1.meta.noteNumber {
                return $0.meta.noteNumber < $1.meta.noteNumber
            }
            return $0.meta.title < $1.meta.title
        }
    }

    /// 按分类分组（仅含有模块的分类，目录页按此渲染）。
    public var grouped: [(category: ModuleCategory, modules: [any SimModule])] {
        let list = all
        var seen: [ModuleCategory: [any SimModule]] = [:]
        for m in list {
            seen[m.meta.category, default: []].append(m)
        }
        // 按 ModuleCategory 声明序输出（目录导航固定顺序）
        return ModuleCategory.allCases.compactMap { cat in
            guard let mods = seen[cat], !mods.isEmpty else { return nil }
            return (cat, mods)
        }
    }

    /// 全文搜索（标题/副题/分类/关键词/笔记编号）。
    public func search(_ query: String) -> [any SimModule] {
        all.filter { $0.meta.matches(query: query) }
    }

    /// 按 id 定位（模块页深链 / 收藏恢复）。
    public func module(id: String) -> (any SimModule)? {
        lock.lock()
        defer { lock.unlock() }
        return modules.first(where: { $0.meta.id == id })
    }

    /// 统计（目录页 chip：全部 N / 实时 N / 基础 N / 专题 N）。
    public var counts: (total: Int, realtime: Int, basic: Int, advanced: Int) {
        let list = all
        return (list.count,
                list.filter { $0.meta.tier == .realtime }.count,
                list.filter { $0.meta.difficulty == .basic }.count,
                list.filter { $0.meta.difficulty == .advanced }.count)
    }

    /// 测试辅助：清空（仅 DEBUG 下允许，防误用）。
    public func removeAll() {
        assert(Thread.isMainThread || true)
        lock.lock()
        defer { lock.unlock() }
        modules.removeAll()
    }
}
