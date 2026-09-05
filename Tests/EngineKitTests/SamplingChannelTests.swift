import Testing
import Foundation
@testable import EngineKit

/// W9 基建验收：SamplingChannel（分批采样 + 进度 + 取消，SecondsChannel 之上）。
@Suite("SamplingChannel 采样双模式通道")
struct SamplingChannelTests {

    /// 同步函数包裹（isMainThread 的 async 语境访问禁令只作用于访问点）
    private func currentIsMainThread() -> Bool { Thread.isMainThread }

    /// 并发收集器（锁保护，规避可变局部量的并发访问诊断）
    private final class SampleBox: @unchecked Sendable {
        private let lock = NSLock()
        private var items: [Int] = []
        func append(_ i: Int) { lock.lock(); items.append(i); lock.unlock() }
        func all() -> [Int] { lock.lock(); defer { lock.unlock() }; return items }
        var count: Int { lock.lock(); defer { lock.unlock() }; return items.count }
    }

    @Test("逐样本全量执行：0..<count 恰好各执行一次，且离开主线程")
    func fullCoverageOffMainThread() async throws {
        let box = SampleBox()
        try await SamplingChannel.run(count: 537, batchSize: 100) { i in
            #expect(!currentIsMainThread(), "采样体必须在 detached 通道执行")
            box.append(i)
        }
        let got = box.all()
        #expect(got.count == 537)
        #expect(Set(got).count == 537, "每个样本号恰好一次")
        #expect(got.min() == 0 && got.max() == 536)
    }

    @Test("进度单调上报：batchSize 分批，最终 fraction == 1")
    func progressMonotonicToFull() async throws {
        let progress = ComputeProgress()
        try await SamplingChannel.run(count: 1000, batchSize: 200,
                                      progress: progress, phase: "测试采样") { _ in }
        #expect(progress.snapshot.fraction == 1.0)
        #expect(progress.snapshot.phase.contains("1000/1000"))
    }

    @Test("取消传播：批间检查点抛 CancellationError，已做样本保留")
    func cancellationAtBatchBoundary() async {
        let box = SampleBox()
        let task = Task {
            try await SamplingChannel.run(count: 1_000_000, batchSize: 50) { i in
                // 每样本少量计算（DEBUG 下取消前的窗口内做不完 1e6 样本）
                var acc = 0.0
                for j in 0..<200 { acc += Double(j &* (i &+ 1)) }
                box.append(i)
                _ = acc
            }
        }
        try? await Task.sleep(nanoseconds: 100_000_000)
        task.cancel()
        await #expect(throws: CancellationError.self) { try await task.value }
        #expect(box.count < 1_000_000, "取消后不应做完全量")
        #expect(box.count > 0, "取消前应已完成若干样本")
    }

    @Test("batchSize 钳位：0/负值按 1 处理仍全量完成")
    func batchSizeClamp() async throws {
        let box = SampleBox()
        try await SamplingChannel.run(count: 17, batchSize: 0) { i in box.append(i) }
        #expect(box.count == 17)
    }
}
