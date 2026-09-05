import Testing
import Foundation
@testable import EngineKit

/// W8 基建验收：秒级调度通道（Task.detached + 取消传播 + 轻量进度）。
@Suite("SecondsChannel 秒级调度通道")
struct SecondsChannelTests {

    /// 同步函数包裹（isMainThread 的 async 语境访问禁令只作用于访问点）
    private func currentIsMainThread() -> Bool { Thread.isMainThread }

    @Test("detached 执行：返回值透传，且不在调用方 actor 上执行（非主线程）")
    func runsDetachedAndReturns() async throws {
        let result = try await SecondsChannel.run { () -> String in
            #expect(!currentIsMainThread(), "重核必须离开调用线程/actor")
            return "ok-\(42)"
        }
        #expect(result == "ok-42")
    }

    @Test("取消传播：外层任务取消 → 作业协作检查点抛 CancellationError")
    func cancellationPropagates() async {
        do {
            let task = Task {
                try await SecondsChannel.run {
                    // 模拟秒级重核：分 50 段，段间协作检查取消
                    var acc = 0.0
                    for i in 0..<50 {
                        try Task.checkCancellation()
                        for j in 0..<2_000_000 { acc += Double(j &* i &+ 1) > 0 ? 1.0 : 0.0 }
                    }
                    return acc
                }
            }
            // 等作业启动后取消（作业整体 > 100ms，此 sleep 落在中段）
            try? await Task.sleep(nanoseconds: 80_000_000)
            task.cancel()
            await #expect(throws: CancellationError.self) { try await task.value }
        }
    }

    @Test("外层已取消时：run 立即取消作业，不上抛结果")
    func alreadyCancelledTaskShortCircuits() async {
        let task = Task {
            try await SecondsChannel.run {
                var acc = 0.0
                for i in 0..<20 {
                    try Task.checkCancellation()
                    for j in 0..<2_000_000 { acc += Double(j &* i &+ 1) > 0 ? 1.0 : 0.0 }
                }
                return acc
            }
        }
        task.cancel()
        await #expect(throws: CancellationError.self) { try await task.value }
    }

    @Test("轻量进度：update 可见 + snapshots 流推最新帧 + finish 收口")
    func progressUpdates() async throws {
        let progress = ComputeProgress()
        let box = SnapshotBox()
        let collector = Task {
            for await snap in progress.snapshots() { box.append(snap) }
        }
        // 模拟 compute 三阶段
        let job = Task.detached {
            for (i, f) in [0.25, 0.5, 1.0].enumerated() {
                progress.update(f, phase: "阶段 \(i + 1)")
                try? await Task.sleep(nanoseconds: 5_000_000)
            }
        }
        _ = await job.value
        try? await Task.sleep(nanoseconds: 20_000_000)
        progress.finish()
        _ = await collector.result

        let received = box.all()
        let final = progress.snapshot
        #expect(final.fraction == 1.0)
        #expect(final.phase == "阶段 3")
        #expect(received.count >= 4, "流应推送初始帧 + 3 次更新，实得 \(received.count)")
        #expect(received.last?.fraction == 1.0)
        // fraction 钳位
        progress.update(7.5)
        #expect(progress.snapshot.fraction == 1.0)
        progress.update(-1)
        #expect(progress.snapshot.fraction == 0.0)
    }

    /// 并发收集器（锁保护的快照盒，规避可变局部量的并发访问诊断）
    private final class SnapshotBox: @unchecked Sendable {
        private let lock = NSLock()
        private var items: [ComputeProgress.Snapshot] = []
        func append(_ s: ComputeProgress.Snapshot) {
            lock.lock(); items.append(s); lock.unlock()
        }
        func all() -> [ComputeProgress.Snapshot] {
            lock.lock(); defer { lock.unlock() }; return items
        }
    }

    @Test("body 抛错：原样上抛，进度流收口（defer finish）")
    func errorPropagation() async {
        struct Probe: Error {}
        let probe = ComputeProgress()
        await #expect(throws: Probe.self) {
            try await SecondsChannel.run(progress: probe) { throw Probe() }
        }
        // finish 后订阅流立即结束（只余已缓冲的最新帧）
        let stream = probe.snapshots()
        var count = 0
        for await _ in stream { count += 1 }
        #expect(count == 1, "finish 后流只剩已缓冲的最新帧")
    }
}
