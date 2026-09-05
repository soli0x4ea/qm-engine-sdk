import Foundation

// MARK: - StrategyCore（W11A 基建）：策略档三件套
//
// 计划 §W11 基建项：降维网格 / 收敛判定 / 缓存命中统计。
// 设计约束：
// - **纯函数/无全局状态**：除 StrategyCache 显式声明为共享缓存外，全部为值语义；
// - **float64 终值**：只处理 Double，CPU 口径（计划 §W11 技术路线）；
// - **可单测**：每个组件独立可用，模块侧组合使用（见 PathIntegralModule）。

// MARK: 降维网格

/// 策略档网格质量档位（计划 §W11：网格 100–200）。
public enum GridQuality: String, Sendable, CaseIterable {
    /// 预览：交互粗算（100 点）
    case preview
    /// 标准：默认口径（150 点）
    case standard
    /// 高精：物理律验证口径（200 点）
    case fine

    /// 该档位的网格点数上限。
    public var gridCount: Int {
        switch self {
        case .preview: return 100
        case .standard: return 150
        case .fine: return 200
        }
    }

    /// 目录/参数面板徽章文案。
    public var badge: String {
        switch self {
        case .preview: return "预览"
        case .standard: return "标准"
        case .fine: return "高精"
        }
    }
}

/// 降维网格策略：参考网格（脚本口径，通常 512–800）→ 策略档计算网格。
public enum StrategyGrid {

    /// 取 quality 档位与参考网格的较小值，并向下取偶（FFT/对称网格友好）。
    ///
    /// 降维前提（SOP §W11）：统计量/收敛量允许比参考网格粗——模块须以
    /// `ConvergenceMonitor` 或物理律测试证明降维后精度仍达标（如路径积分
    /// 200 点已达到 Mehler relerr 5.1e-4 的方案验证口径）。
    public static func reduced(reference: Int, quality: GridQuality) -> Int {
        let n = min(max(reference, 2), quality.gridCount)
        return n % 2 == 0 ? n : n - 1
    }

    /// 等距抽稀到目标点数（保端点；纯函数，用于把计算网格曲线降到显示分辨率）。
    public static func downsample(_ ys: [Double], to target: Int) -> [Double] {
        precondition(target >= 2, "downsample: target >= 2")
        guard ys.count > target else { return ys }
        var out = [Double]()
        out.reserveCapacity(target)
        let last = ys.count - 1
        for k in 0..<target {
            let i = k == target - 1 ? last : (k * last) / (target - 1)
            out.append(ys[i])
        }
        return out
    }
}

// MARK: 收敛判定

/// 标量序列收敛监视器：连续 `patience` 步相对变化 < `relTolerance` 即判收敛；
/// 步数达到 `maxSteps` 强制停机（返回 shouldStop = true，但不标记 convergedAt）。
///
/// 典型用法（策略档迭代/加密循环）：
/// ```swift
/// var monitor = ConvergenceMonitor(relTolerance: 2e-3, patience: 2, maxSteps: 10)
/// for step in 0..<10 {
///     let value = estimate(step)
///     if monitor.observe(value) { break }   // 收敛或到步数上限
/// }
/// ```
public struct ConvergenceMonitor: Sendable {
    public let relTolerance: Double
    public let patience: Int
    public let maxSteps: Int
    /// 已观测步数。
    public private(set) var steps = 0
    /// 容差收敛发生的步号（1 起；nil = 尚未收敛）。
    public private(set) var convergedAt: Int?
    private var lastValue: Double?
    private var stableRun = 0

    public init(relTolerance: Double, patience: Int = 2, maxSteps: Int) {
        precondition(relTolerance > 0, "relTolerance 必须为正")
        precondition(patience >= 1, "patience >= 1")
        precondition(maxSteps >= 1, "maxSteps >= 1")
        self.relTolerance = relTolerance
        self.patience = patience
        self.maxSteps = maxSteps
    }

    /// 是否已容差收敛。
    public var isConverged: Bool { convergedAt != nil }

    /// 观测新值：返回 true 表示应停机（容差收敛或达到 maxSteps）。
    public mutating func observe(_ value: Double) -> Bool {
        steps += 1
        defer { lastValue = value }
        if steps >= maxSteps {
            // 最后一步：先判收敛再报停机
            return check(value) || true
        }
        return check(value)
    }

    private mutating func check(_ value: Double) -> Bool {
        guard let last = lastValue, last != 0 else { return false }
        let rel = abs(value - last) / abs(last)
        stableRun = rel < relTolerance ? stableRun + 1 : 0
        if stableRun >= patience {
            if convergedAt == nil { convergedAt = steps }
            return true
        }
        return false
    }
}

// MARK: 参数哈希缓存（命中统计 + LRU）

/// 缓存统计快照（线程安全读取）。
public struct StrategyCacheStats: Sendable, Equatable {
    public private(set) var hits: Int
    public private(set) var misses: Int
    public private(set) var evictions: Int

    public var lookups: Int { hits + misses }
    /// 命中率（未发生查找时为 0）。
    public var hitRate: Double { lookups == 0 ? 0 : Double(hits) / Double(lookups) }

    init(hits: Int = 0, misses: Int = 0, evictions: Int = 0) {
        self.hits = hits
        self.misses = misses
        self.evictions = evictions
    }

    /// 内部计数（StrategyCache 持锁调用）。
    mutating func recordHit() { hits += 1 }
    mutating func recordMiss() { misses += 1 }
    mutating func recordEviction() { evictions += 1 }
}

/// 参数哈希缓存：策略档重型中间结果（如特征分解）按参数键复用，
/// 带**命中/未命中/逐出**统计——周报「缓存命中统计」与模块摘要的数据源。
///
/// - 线程安全（NSLock）；LRU 逐出，`capacity` 上界防内存漂移；
/// - 与 NSCache 的差异：统计可测、逐出确定（单测可断言），故自持小字典而非 NSCache
///   （计划 §3.1 的 NSCache 参数哈希缓存模式在本层形式化）。
public final class StrategyCache<Key: Hashable & Sendable, Value: Sendable>: @unchecked Sendable {

    private let lock = NSLock()
    private var storage: [Key: Value] = [:]
    private var order: [Key] = []          // 队首 = 最近使用
    private let capacity: Int
    private var statsValue = StrategyCacheStats()

    public init(capacity: Int = 8) {
        precondition(capacity >= 1, "capacity >= 1")
        self.capacity = capacity
    }

    /// 取缓存值：命中计 hits，未命中计 misses。
    public func value(for key: Key) -> Value? {
        lock.lock()
        defer { lock.unlock() }
        if let v = storage[key] {
            statsValue.recordHit()
            touch(key)
            return v
        }
        statsValue.recordMiss()
        return nil
    }

    /// 插入/覆盖；容量超限时 LRU 逐出并计 evictions。
    public func insert(_ value: Value, for key: Key) {
        lock.lock()
        defer { lock.unlock() }
        if storage[key] == nil, storage.count >= capacity {
            let evicted = order.removeLast()
            storage[evicted] = nil
            statsValue.recordEviction()
        }
        storage[key] = value
        touch(key)
    }

    /// 当前统计快照。
    public var stats: StrategyCacheStats {
        lock.withLock { statsValue }
    }

    /// 清空统计（保留数据）。
    public func resetStats() {
        lock.withLock { statsValue = StrategyCacheStats() }
    }

    /// 清空数据与统计。
    public func removeAll() {
        lock.withLock {
            storage.removeAll()
            order.removeAll()
            statsValue = StrategyCacheStats()
        }
    }

    /// LRU 触碰：移到队首。
    private func touch(_ key: Key) {
        order.removeAll { $0 == key }
        order.insert(key, at: 0)
    }
}

// MARK: 参数哈希键

/// 稳定参数哈希键：模块 id + 有序 (名称, Double) 对。
/// 浮点按位十六进制编码（bitPattern），跨进程/跨调用位级稳定，无格式化歧义。
public enum StrategyKey {

    public static func make(_ moduleID: String, _ params: [(String, Double)]) -> String {
        let body = params
            .map { "\($0.0)=\(String(format: "%016llx", $0.1.bitPattern))" }
            .joined(separator: ";")
        return "\(moduleID)#\(body)"
    }
}
