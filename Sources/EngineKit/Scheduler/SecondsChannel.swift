import Foundation

// MARK: - 秒级调度通道（W8 基建）

/// 秒级计算进度（线程安全，轻量）。
///
/// 模块 compute 内部按阶段 `update(_:phase:)`；UI 侧两种消费方式（W8 交付通道能力，
/// App 面板接线随 iPad/性能周落地）：
/// - 轮询 `snapshot`（如 Combine/Timer 采样）；
/// - 订阅 `snapshots()`（AsyncStream，缓冲最新一帧，背压安全）。
public final class ComputeProgress: @unchecked Sendable {

    /// 进度快照：fraction ∈ [0, 1]，phase 为阶段短语（可空串）。
    public struct Snapshot: Sendable, Equatable {
        public let fraction: Double
        public let phase: String
    }

    private let lock = NSLock()
    private var current = Snapshot(fraction: 0, phase: "")
    private var continuations: [Int: AsyncStream<Snapshot>.Continuation] = [:]
    private var nextID = 0
    private var finished = false

    public init() {}

    /// 更新进度（fraction 自动钳到 [0,1]；phase 为 nil 时保留原阶段文案）。
    public func update(_ fraction: Double, phase: String? = nil) {
        let snap = Snapshot(fraction: min(max(fraction, 0), 1),
                            phase: phase ?? lock.withLock { current.phase })
        var sinks: [AsyncStream<Snapshot>.Continuation] = []
        lock.lock()
        current = snap
        sinks = Array(continuations.values)
        lock.unlock()
        for c in sinks { c.yield(snap) }
    }

    /// 当前快照。
    public var snapshot: Snapshot {
        lock.withLock { current }
    }

    /// 订阅进度流（compute 收尾时 SecondsChannel 自动 `finish()`；晚到订阅立即收口）。
    public func snapshots() -> AsyncStream<Snapshot> {
        AsyncStream(bufferingPolicy: .bufferingNewest(1)) { continuation in
            lock.lock()
            let isFinished = finished
            continuation.yield(current)
            if isFinished {
                continuation.finish()
            } else {
                continuations[nextID] = continuation
                nextID += 1
            }
            lock.unlock()
        }
    }

    /// 结束全部订阅（compute 正常/异常收尾时调用；SecondsChannel 已代为调用）。
    public func finish() {
        lock.lock()
        finished = true
        let sinks = Array(continuations.values)
        continuations.removeAll()
        lock.unlock()
        for c in sinks { c.finish() }
    }
}

/// 秒级档调度通道：`Task.detached` + 取消传播 + 轻量进度。
///
/// 设计（计划 §W8 基建项 / §3.2 秒级档约束）：
/// - **detached**：模块 compute 常由 UI 侧 Task（MainActor 语境）调用——重核
///   （数千阶 eigh、万点扫描）必须搬离调用方 actor，主线程只等结果，交互不掉帧；
/// - **取消传播**：外层 Task 取消 → `onCancel` 立即取消 detached 作业；
///   作业体内按 SOP 周期性 `try Task.checkCancellation()`（协作检查点），
///   取消以 `CancellationError` 上抛，ModuleRunner 按「被新计算取代」静默丢弃；
/// - **轻量进度**：`ComputeProgress` 锁保护 + AsyncStream 最新帧，开销 ~纳秒级，
///   不进 compute 热路径的分配器。
///
/// 用法（tier-3 模块 compute 内）：
/// ```swift
/// let progress = ComputeProgress()
/// return try await SecondsChannel.run(progress: progress) {
///     progress.update(0.1, phase: "组装哈密顿量")
///     ... // 阶段间 try Task.checkCancellation()
/// }
/// ```
public enum SecondsChannel {

    /// 在优先级 `priority` 的后台任务中执行 `body`；外层任务取消时作业连带取消。
    ///
    /// - Parameters:
    ///   - progress: 可选进度句柄（body 内捕获使用）；作业收尾自动 `finish()`。
    ///   - body: 计算核。Sendable 闭包，纯函数、无全局状态（SOP 约束）。
    public static func run<T: Sendable>(
        priority: TaskPriority = .userInitiated,
        progress: ComputeProgress? = nil,
        _ body: @escaping @Sendable () async throws -> T
    ) async throws -> T {
        let job = Task.detached(priority: priority) {
            try await body()
        }
        defer { progress?.finish() }
        return try await withTaskCancellationHandler {
            try await job.value
        } onCancel: {
            job.cancel()
        }
    }
}
