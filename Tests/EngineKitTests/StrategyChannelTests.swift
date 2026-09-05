import Testing
import Foundation
@testable import EngineKit

/// 策略档调度通道单测（W11A 基建）：串行互斥 / 进度与取消 / 检查点续算。
@Suite(.serialized)
struct StrategyChannelTests {

    /// 作业时间线记录（线程安全）。
    private actor Timeline {
        var entries: [(id: String, at: ContinuousClock.Instant)] = []
        func mark(_ id: String) { entries.append((id, ContinuousClock.now)) }
        func intervals() -> [(String, ContinuousClock.Instant, ContinuousClock.Instant)] {
            var out: [(String, ContinuousClock.Instant, ContinuousClock.Instant)] = []
            var open: (String, ContinuousClock.Instant)?
            for e in entries {
                if e.id.hasSuffix(":enter") {
                    open = (String(e.id.dropLast(6)), e.at)
                } else if e.id.hasSuffix(":leave"), let o = open {
                    out.append((o.0, o.1, e.at))
                    open = nil
                }
            }
            return out
        }
    }

    @Test("串行互斥：并发两个作业时间线不重叠（FIFO 过门）")
    func serialExecution() async throws {
        let timeline = Timeline()
        func job(_ id: String, ms: UInt64) async throws {
            _ = try await StrategyChannel.run(session: StrategySession<Int>(), progress: nil) { (_: StrategyContext<Int>) in
                await timeline.mark("\(id):enter")
                try? await Task.sleep(nanoseconds: ms * 1_000_000)
                await timeline.mark("\(id):leave")
            }
        }
        async let a: Void = job("A", ms: 120)
        async let b: Void = job("B", ms: 20)
        _ = try await (a, b)
        let intervals = await timeline.intervals()
        #expect(intervals.count == 2)
        // 区间不重叠：后者 enter ≥ 前者 leave
        let sorted = intervals.sorted { $0.1 < $1.1 }
        #expect(sorted[1].1 >= sorted[0].2, "B 必须等 A 释放门后才 enter")
    }

    @Test("取消传播：外层取消 → 作业以 CancellationError 收场")
    func cancellation() async throws {
        let progress = ComputeProgress()
        let task = Task {
            try await StrategyChannel.run(session: StrategySession<Int>(), progress: progress) { (_: StrategyContext<Int>) in
                progress.update(0.5, phase: "running")
                while !Task.isCancelled {
                    try await Task.sleep(nanoseconds: 5_000_000)
                }
                try Task.checkCancellation()
                return 0
            }
        }
        // 等作业启动（进度 ≥ 0.5）再取消
        while progress.snapshot.fraction < 0.5 {
            try await Task.sleep(nanoseconds: 2_000_000)
        }
        task.cancel()
        await #expect(throws: CancellationError.self) { try await task.value }
        #expect(progress.snapshot.fraction >= 0.5)
    }

    @Test("检查点续算：取消后同会话恢复，resumedFrom = 上次检查点")
    func checkpointResume() async throws {
        let session = StrategySession<Int>()
        // 第一次 run：落检查点后人为抛取消（模拟作业中断）
        await #expect(throws: CancellationError.self) {
            try await StrategyChannel.run(session: session) { ctx in
                ctx.checkpoint(2)
                throw CancellationError()
            }
        }
        #expect(session.current == 2, "取消后检查点保留在会话")
        // 第二次 run：从检查点恢复
        let resumed = try await StrategyChannel.run(session: session) { ctx -> Int in
            ctx.checkpoint(4)
            return ctx.resumedFrom ?? -1
        }
        #expect(resumed == 2, "resumedFrom = 上次会话检查点")
        #expect(session.current == 4)
    }

    @Test("排队期取消：门被占时外层取消 → 不启动作业且门正常释放")
    func cancelWhileQueued() async throws {
        let session = StrategySession<Int>()
        // 占门的长作业
        let blocker = Task {
            try await StrategyChannel.run(session: session) { _ in
                try await Task.sleep(nanoseconds: 200_000_000)
                return 0
            }
        }
        try await Task.sleep(nanoseconds: 20_000_000)   // 确保 blocker 已持门
        // 排队者随即取消
        let queued = Task {
            try await StrategyChannel.run(session: session) { _ in
                // 无论取消发生在排队期还是作业体，sleep 都会以取消收场
                try await Task.sleep(nanoseconds: 500_000_000)
                return 1
            }
        }
        queued.cancel()
        await #expect(throws: CancellationError.self) { try await queued.value }
        _ = try await blocker.value        // 门正常释放
        // 门回收后新作业可正常进（含恢复点）
        let ok = try await StrategyChannel.run(session: session) { ctx -> Int in
            ctx.checkpoint(1)
            return 42
        }
        #expect(ok == 42)
    }
}
