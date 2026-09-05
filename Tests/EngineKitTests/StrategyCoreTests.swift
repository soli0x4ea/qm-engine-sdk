import Testing
import Foundation
@testable import EngineKit

/// StrategyCore 三件套单测（W11A 基建）：降维网格 / 收敛判定 / 参数哈希缓存。
@Suite(.serialized)
struct StrategyCoreTests {

    // MARK: 降维网格

    @Test("StrategyGrid：档位上限 + 参考网格取小 + 保偶数")
    func gridReduction() {
        #expect(StrategyGrid.reduced(reference: 800, quality: .preview) == 100)
        #expect(StrategyGrid.reduced(reference: 800, quality: .standard) == 150)
        #expect(StrategyGrid.reduced(reference: 800, quality: .fine) == 200)
        // 参考网格比档位还小：不放大
        #expect(StrategyGrid.reduced(reference: 120, quality: .standard) == 120)
        #expect(StrategyGrid.reduced(reference: 64, quality: .fine) == 64)
        // 保偶（奇数档位向下取偶）
        #expect(StrategyGrid.reduced(reference: 800, quality: .preview) % 2 == 0)
        // 退化输入
        #expect(StrategyGrid.reduced(reference: 1, quality: .fine) >= 2)
    }

    @Test("StrategyGrid.downsample：保端点、目标点数、超短不动")
    func downsample() {
        let ys = Array(0..<101).map(Double.init)
        let out = StrategyGrid.downsample(ys, to: 51)
        #expect(out.count == 51)
        #expect(out.first == 0 && out.last == 100)
        // 单调不减（保序抽稀）
        #expect(zip(out, out.dropFirst()).allSatisfy { $0 <= $1 })
        // 目标 ≥ 原点数：原样返回
        #expect(StrategyGrid.downsample(ys, to: 200) == ys)
    }

    // MARK: 收敛判定

    @Test("ConvergenceMonitor：几何收敛序列在 patience 步后停机")
    func convergenceOnGeometric() {
        var m = ConvergenceMonitor(relTolerance: 1e-2, patience: 2, maxSteps: 50)
        // E₀ 型收敛序列：0.73 → 0.538 → 0.507 → 0.501 → 0.500 …
        let seq = [0.7309, 0.5378, 0.5072, 0.5012, 0.5000, 0.4999, 0.4998]
        var stopped = false
        for v in seq where !stopped {
            stopped = m.observe(v)
        }
        #expect(m.isConverged)
        #expect(m.convergedAt != nil)
        #expect(m.steps <= seq.count)
    }

    @Test("ConvergenceMonitor：噪声序列不误判，maxSteps 强制停机")
    func maxStepsCap() {
        var m = ConvergenceMonitor(relTolerance: 1e-6, patience: 2, maxSteps: 5)
        var stoppedAt = -1
        for k in 1...20 {
            let v = 0.5 + 0.1 * sin(Double(k))
            if m.observe(v) { stoppedAt = k; break }
        }
        #expect(!m.isConverged, "振荡序列不满足容差收敛")
        #expect(stoppedAt == 5, "第 5 步（maxSteps）强制停机")
    }

    @Test("ConvergenceMonitor：patience = 1 单步即判")
    func patienceOne() {
        var m = ConvergenceMonitor(relTolerance: 0.1, patience: 1, maxSteps: 10)
        let first = m.observe(1.0)
        #expect(!first, "首值无前值可比，不判")
        let second = m.observe(1.05)
        #expect(second, "单步 rel 0.05 < 0.1 即收敛")
        #expect(m.convergedAt == 2)
    }

    // MARK: 参数哈希缓存

    @Test("StrategyCache：miss→insert→hit 统计与 LRU 逐出")
    func cacheStatsAndLRU() {
        let cache = StrategyCache<String, Int>(capacity: 2)
        #expect(cache.value(for: "a") == nil)      // miss
        cache.insert(1, for: "a")
        cache.insert(2, for: "b")
        cache.insert(3, for: "c")                  // 逐出 LRU "a"
        #expect(cache.value(for: "a") == nil)      // miss（被逐出）
        #expect(cache.value(for: "b") == 2)        // hit
        #expect(cache.value(for: "c") == 3)        // hit
        let stats = cache.stats
        #expect(stats.hits == 2)
        #expect(stats.misses == 2)
        #expect(stats.evictions == 1)
        #expect(abs(stats.hitRate - 0.5) < 1e-12)
        cache.resetStats()
        #expect(cache.stats == StrategyCacheStats())
        #expect(cache.value(for: "b") == 2, "resetStats 不清数据")
        cache.removeAll()
        #expect(cache.value(for: "b") == nil)
    }

    @Test("StrategyCache：覆盖不逐出、容量上界不漂移")
    func cacheOverwrite() {
        let cache = StrategyCache<String, Int>(capacity: 2)
        cache.insert(1, for: "a")
        cache.insert(11, for: "a")                 // 覆盖，不逐出
        #expect(cache.value(for: "a") == 11)
        #expect(cache.stats.evictions == 0)
        #expect(cache.stats.hits == 1, "仅 value(for:) 计命中")
        cache.insert(2, for: "b")
        cache.insert(3, for: "c")
        #expect(cache.stats.evictions == 1, "容量 2，第三个键逐出 LRU（此时 LRU 为 b）")
    }

    // MARK: 参数哈希键

    @Test("StrategyKey：同参同键、异参异键、位级稳定")
    func keyStability() {
        let k1 = StrategyKey.make("模块A", [("grid", 200), ("T", 4.0)])
        let k2 = StrategyKey.make("模块A", [("grid", 200), ("T", 4.0)])
        #expect(k1 == k2)
        #expect(k1 != StrategyKey.make("模块A", [("grid", 200), ("T", 5.0)]))
        #expect(k1 != StrategyKey.make("模块A", [("grid", 150), ("T", 4.0)]))
        #expect(k1 != StrategyKey.make("模块B", [("grid", 200), ("T", 4.0)]))
        // 顺序敏感（约定固定序）
        #expect(StrategyKey.make("M", [("a", 1), ("b", 2)])
                != StrategyKey.make("M", [("b", 2), ("a", 1)]))
        // 浮点位级稳定（0.1 的二进制表示跨调用一致）
        #expect(StrategyKey.make("M", [("x", 0.1)])
                == StrategyKey.make("M", [("x", 0.1)]))
    }
}
