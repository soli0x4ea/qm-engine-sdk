import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// 模块 fixtures 统一加载器：fixtures/modules/<脚本>__v2022.json（W1 管线产物）。
/// 结构：figures[] → axes[] → { lines[{label,x,y}], bars[{x,height,width}] }；asserts[{line,expr,ok}]。
struct ModuleFixture {
    let document: [String: Any]

    static func load(_ script: String) throws -> ModuleFixture {
        let base = (script as NSString).deletingPathExtension
        let url = try #require(
            Bundle.module.url(forResource: base, withExtension: "json", subdirectory: "Fixtures"))
        let data = try Data(contentsOf: url)
        return ModuleFixture(document: try #require(
            try JSONSerialization.jsonObject(with: data) as? [String: Any]))
    }

    private func axis(_ figure: Int, _ axisIndex: Int) throws -> [String: Any] {
        let figures = try #require(document["figures"] as? [[String: Any]])
        let axes = try #require(figure < figures.count ? figures[figure]["axes"] as? [[String: Any]] : nil)
        return try #require(axisIndex < axes.count ? axes[axisIndex] : nil)
    }

    /// 取指定标签曲线的 (x, y)。
    func line(_ figure: Int, _ axisIndex: Int, label: String) throws -> (x: [Double], y: [Double]) {
        let ax = try axis(figure, axisIndex)
        let lines = try #require(ax["lines"] as? [[String: Any]])
        let ln = try #require(lines.first { ($0["label"] as? String)?.contains(label) == true })
        return (try #require((ln["x"] as? [NSNumber])?.map(\.doubleValue)),
                try #require((ln["y"] as? [NSNumber])?.map(\.doubleValue)))
    }

    /// 按下标取曲线（标签子串歧义时用，如 "nonrelativistic" 包含 "relativistic"）。
    func line(_ figure: Int, _ axisIndex: Int, index: Int) throws -> (x: [Double], y: [Double]) {
        let ax = try axis(figure, axisIndex)
        let lines = try #require(ax["lines"] as? [[String: Any]])
        let ln = try #require(index < lines.count ? lines[index] : nil)
        return (try #require((ln["x"] as? [NSNumber])?.map(\.doubleValue)),
                try #require((ln["y"] as? [NSNumber])?.map(\.doubleValue)))
    }

    /// 散点数据（无标签，按出现序）：[(x, y)] 每组一条。
    func scatters(_ figure: Int, _ axisIndex: Int) throws -> [(x: [Double], y: [Double])] {
        let ax = try axis(figure, axisIndex)
        let raw = ax["scatters"] as? [[String: Any]] ?? []
        return try raw.map {
            (try #require(($0["x"] as? [NSNumber])?.map(\.doubleValue)),
             try #require(($0["y"] as? [NSNumber])?.map(\.doubleValue)))
        }
    }

    /// 柱状数据：[(x, height, width)]。
    func bars(_ figure: Int, _ axisIndex: Int) throws -> [(x: Double, height: Double, width: Double)] {
        let ax = try axis(figure, axisIndex)
        let raw = try #require(ax["bars"] as? [[String: Any]])
        return try raw.map {
            (try #require(($0["x"] as? NSNumber)?.doubleValue),
             ($0["height"] as? NSNumber)?.doubleValue ?? 0,
             ($0["width"] as? NSNumber)?.doubleValue ?? 0)
        }
    }
}

/// 逐点相对误差断言（fixtures 存 9 位有效数字，阈值取 1e-8）。
func expectPointwiseClose(
    _ actual: [Double], _ expected: [Double], tolerance: Double,
    _ context: @autoclosure () -> String, sourceLocation: SourceLocation = #_sourceLocation
) {
    #expect(actual.count == expected.count, "\(context()): 点数 \(actual.count) vs \(expected.count)",
            sourceLocation: sourceLocation)
    guard actual.count == expected.count else { return }
    for i in actual.indices {
        let denom = max(abs(expected[i]), 1e-300)
        #expect(abs(actual[i] - expected[i]) / denom < tolerance,
                "\(context()) 第 \(i) 点：\(actual[i]) vs \(expected[i])",
                sourceLocation: sourceLocation)
    }
}

/// 性能断言：耗时 < 预算毫秒（实时档 16 ms / 秒级档 2000 ms）。
///
/// **计时口径：best-of-N 最小值（墙钟）。**
///
/// Swift Testing 默认并行执行套件，单次墙钟采样会把其余测试线程与系统负载的
/// 等待时间计入 compute 成本。实测同一模块（W4 指针耦合与芝诺）串行执行 3.96 ms，
/// 全量并行下中位被放大到 25–31 ms——即噪声可达真值的 7 倍。
///
/// 争用只会**抬高**耗时、不会压低，故 N 次采样的最小值是无争用成本的一致估计
/// （微基准在共享机器上的通行做法 best-of-N）；中位/均值则随并行度漂移。
/// 曾尝试 `CLOCK_THREAD_CPUTIME_ID` 线程 CPU 时间，但 Darwin 实现粒度粗、
/// 同样被调度放大（实测 28 ms），故不采用。
///
/// 为压低采样噪声：先 `warmup` 次预热（稳定 JIT/缓存/对象池），再取 `runs` 次
/// 样本的最小值作为判据。样本越多，最小值越逼近真实无争用成本。
/// 中位一并打印，便于人工核对真实交互延迟与争用程度。
///
/// **超时不阻断政策（2026-09-05 起）**：`enforce = false` 时只测量并打印，不做断言——
/// 超预算模块登记进周报告的「超预算登记表」，统一留到全部开发完成后的性能优化专项处理，
/// 不再单周就地优化（避免为压耗时牺牲数值口径）。既有 41 个模块保持 `enforce = true`。
///
/// **收尾 1b 提升（2026-09-05）**：全量 63 模块性能审计完成——零超预算、零 >50% 预算
/// （实时档最重 3.91 ms / 16，秒级档最重 104.9 ms / 2000，策略档最重 22.8 ms / 10000），
/// 原 16 处 `enforce = false` 欠账测试全部提升为强制口径（enforce = true）。
/// 全量账本见 build/收尾执行包1b_全量性能审计/ 周报告。
func expectComputeUnderBudget(
    module: any SimModule, values: ParamValues, constants: ConstantsSet,
    budgetMillis: Double, runs: Int = 41, warmup: Int = 5, enforce: Bool = true,
    sourceLocation: SourceLocation = #_sourceLocation
) async throws {
    for _ in 0..<warmup {
        _ = try await module.compute(values, constants: constants)
    }
    var samples: [Double] = []
    for _ in 0..<runs {
        let t0 = ContinuousClock.now
        _ = try await module.compute(values, constants: constants)
        let d = ContinuousClock.now - t0
        let c = d.components
        samples.append(Double(c.seconds) * 1000 + Double(c.attoseconds) / 1e15)
    }
    samples.sort()
    let best = samples[0]
    let median = samples[samples.count / 2]
    let overBudget = best >= budgetMillis
    let tag = enforce ? (overBudget ? "❌超预算" : "✅") : (overBudget ? "⚠️超预算·仅登记" : "✅")
    print("⏱ \(module.meta.id) compute 最小 \(String(format: "%.2f", best)) ms"
          + "（best-of-\(runs)，预算 \(budgetMillis) ms）· 中位 \(String(format: "%.2f", median)) ms · \(tag)")
    guard enforce else { return }
    let detailText = "compute 耗时 \(String(format: "%.2f", best)) ms 超预算 \(budgetMillis) ms"
        + "（中位 \(String(format: "%.2f", median)) ms 含并行争用，判据取 best-of-\(runs) 最小值）"
    #expect(best < budgetMillis, Comment(stringLiteral: detailText),
            sourceLocation: sourceLocation)
}
