import Foundation

// MARK: - SamplingChannel（W9 基建：采样双模式）

/// 采样双模式契约（计划 §W9：交互采样数 vs 完整后台完成）。
///
/// - `interactive`：滑杆指定的小样本量（如纠缠模块 200–500），交互拖动即时反馈；
/// - `full`：Python 脚本同规模的全量样本（如 4000），经 `SamplingChannel.run`
///   脱离主线程分批执行，带进度与取消，UI 保持响应。
///
/// 种子口径见 docs/RNG_SEED_POLICY.md：种子 = Python 脚本原值（不作为用户参数暴露），
/// Swift 侧 SplitMix64 + Box-Muller 同分布复采，统计断言对拍。
public enum SampleMode: String, Sendable {
    /// 交互模式：样本数由滑杆给出（快速反馈）
    case interactive
    /// 完整后台：固定全量样本数（脚本口径），分批 + 进度 + 可取消
    case full
}

/// 采样任务的执行通道：在 `SecondsChannel` 之上加「分批循环 + 进度上报 + 取消检查」。
///
/// 用法（秒级档模块 compute 内）：
///
///     let result = try await SamplingChannel.run(
///         count: 4000, batchSize: 200, progress: progress
///     ) { index, ctx in
///         ctx.accumulate(sample(index))   // 或直接在 body 里逐样本处理
///     }
///
/// 语义：
/// - 整个 body 在 `Task.detached` 上执行（复用 SecondsChannel，离开主线程）；
/// - 每 `batchSize` 个样本调用一次 `progress.update(_, phase:)` 并检查取消；
/// - `ComputeProgress.finish()` 由 SecondsChannel 的 defer 兜底。
public enum SamplingChannel {

    /// 分批采样执行器。
    ///
    /// - Parameters:
    ///   - count: 总样本数
    ///   - batchSize: 进度上报 / 取消检查的批大小（默认 200）
    ///   - priority: detached 任务优先级（默认 .userInitiated）
    ///   - progress: 进度对象（nil = 不上报，仍检查取消）
    ///   - phase: 进度文案前缀（如 "Ginibre 采样"）
    ///   - body: 单样本处理；入参为样本序号（0..<count）。抛错即整体失败。
    public static func run(
        count: Int,
        batchSize: Int = 200,
        priority: TaskPriority = .userInitiated,
        progress: ComputeProgress? = nil,
        phase: String = "采样",
        _ body: @escaping @Sendable (Int) throws -> Void
    ) async throws -> Void {
        let batch = max(1, batchSize)
        try await SecondsChannel.run(priority: priority, progress: progress) {
            for i in 0..<count {
                try body(i)
                if i % batch == batch - 1 || i == count - 1 {
                    try Task.checkCancellation()
                    progress?.update(Double(i + 1) / Double(count),
                                     phase: "\(phase) \(i + 1)/\(count)")
                }
            }
        }
    }
}
