import Foundation

// MARK: - 策略档调度通道（W11A 基建）
//
// 计划 §W11：调度器策略通道（串行队列 + 进度/取消 + 恢复进度）。
// 与 SecondsChannel/SamplingChannel 的分工：
// - **串行互斥**：全局 `StrategyGate` 保证同一时刻只有一个策略档作业在跑，
//   重型计算不挤占实时档交互（计划 §W11 两带光学模型「后台队列串行防挤占」）；
// - **检查点续算**：作业体在阶段边界 `ctx.checkpoint(_:)` 落盘进度状态，
//   被取消后再次 `run` 同一 `StrategySession` 时从 `ctx.resumedFrom` 恢复，
//   已完成阶段可跳过——对应「后台续算/恢复进度」的最小内核。
//   （BGTaskScheduler 注册属 iOS App 层，按 W11A 范围不在引擎内实现。）
// - **优先级**：`.utility`（策略档后台语义，低于秒级档 `.userInitiated`）。

/// 策略档串行门：同一时刻至多持有一个作业（actor 语义保证）。
public actor StrategyGate {

    /// 全局共享门：所有策略档作业经此串行。
    public static let shared = StrategyGate()

    private var locked = false
    private var waiters: [CheckedContinuation<Void, Never>] = []

    /// 获取门（已被持有则挂起排队；FIFO 唤醒）。
    public func acquire() async {
        if locked {
            await withCheckedContinuation { waiters.append($0) }
            // 唤醒即持有：release 移交所有权时不复位 locked
        } else {
            locked = true
        }
    }

    /// 释放门：有等待者则移交所有权（保持 locked），否则复位。
    public func release() {
        if let next = waiters.first {
            waiters.removeFirst()
            next.resume()
        } else {
            locked = false
        }
    }
}

/// 可续算会话：持有跨 run 的检查点状态（线程安全）。
///
/// 泛型 `State` 为作业自定义的轻量进度标记（如已完成阶段号、网格游标），
/// 须 Sendable；重型数据本体走 StrategyCache，不塞进会话。
public final class StrategySession<State: Sendable>: @unchecked Sendable {

    private let lock = NSLock()
    private var state: State?

    public init(initial: State? = nil) {
        self.state = initial
    }

    /// 当前检查点（nil = 尚无）。
    public var current: State? {
        lock.withLock { state }
    }

    /// 记录检查点（作业体内经 StrategyContext.checkpoint 调用）。
    public func setCheckpoint(_ newState: State) {
        lock.withLock { state = newState }
    }

    /// 重置会话（测试用）。
    public func reset(to newState: State? = nil) {
        lock.withLock { state = newState }
    }
}

/// 传给作业体的上下文：恢复点 + 检查点写入。
public struct StrategyContext<State: Sendable>: Sendable {

    /// 本次 run 的恢复起点（会话此前记录的检查点；首次为 nil）。
    public let resumedFrom: State?
    private let session: StrategySession<State>

    init(resumedFrom: State?, session: StrategySession<State>) {
        self.resumedFrom = resumedFrom
        self.session = session
    }

    /// 落盘检查点（取消后仍保留在会话中，供下次 run 恢复）。
    public func checkpoint(_ state: State) {
        session.setCheckpoint(state)
    }
}

/// 策略档调度通道：串行门 + detached 作业 + 取消传播 + 轻量进度 + 检查点续算。
public enum StrategyChannel {

    /// 全局串行门（策略档作业互斥）。
    public static let gate = StrategyGate.shared

    /// 在策略档通道中执行 `body`。
    ///
    /// - Parameters:
    ///   - session: 可续算会话（nil 时用一次性会话，不跨 run 保留检查点）。
    ///   - progress: 可选进度句柄（body 内捕获使用；收尾自动 finish）。
    ///   - body: 计算核。Sendable 闭包，收到 `StrategyContext`；
    ///     体内应周期性 `try Task.checkCancellation()` 协作取消。
    ///
    /// 取消语义与 SecondsChannel 一致：外层 Task 取消 → 作业连带取消，
    /// 以 `CancellationError` 上抛；会话检查点已落盘的部分不丢失。
    public static func run<State: Sendable, T: Sendable>(
        session: StrategySession<State>? = nil,
        progress: ComputeProgress? = nil,
        _ body: @escaping @Sendable (StrategyContext<State>) async throws -> T
    ) async throws -> T {
        let own = session ?? StrategySession<State>()
        let ctx = StrategyContext<State>(resumedFrom: own.current, session: own)

        await gate.acquire()
        defer { progress?.finish() }
        do {
            try Task.checkCancellation()   // 排队期间被取消：不启动作业
            let job = Task.detached(priority: .utility) {
                try await body(ctx)
            }
            let value = try await withTaskCancellationHandler {
                try await job.value
            } onCancel: {
                job.cancel()
            }
            await gate.release()
            return value
        } catch {
            await gate.release()
            throw error
        }
    }
}
